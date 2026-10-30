from __future__ import annotations
import argparse
import pandas as pd
from pathlib import Path
import pysam
import json
import math
from dataclasses import dataclass,field
SVTYPE_SUFFIX={"INS","DEL"}
CHROM_TAG=["SNP","INS","DEL","ALL"]
TOLERANCE_RATE=[0,0.2,0.4,0.6,0.8,1]
@dataclass
class Truth:
    id: int
    pos: int
    start:int | None= None
    end: int | None= None
    matched_candidates: list[str] = field(default_factory=list)

@dataclass
class Candidate:
    id: str
    start: int
    end: int
    
truth_groups = {
    tag: [] for tag in CHROM_TAG
}

candidate_groups = {
    tag: [] for tag in CHROM_TAG
}

result={
    #"candidate_total":0,
    #"candidate_tp":0,
    #"candidate_fp":0,
    #"truth_total":0,
    #"truth_tp":0,
    #"truth_fn":0,
    #"precision":0,
    #"recall":0,
    #"F1":0
}
ROOT = Path(__file__).resolve().parent
if ROOT.name == "analysis":
    ROOT = ROOT.parent

def project_path(path: str) -> Path:
    """Interpret relative paths from the project root, not the current cwd."""
    p = Path(path).expanduser()
    return (p if p.is_absolute() else ROOT / p).resolve()

def candidate_count(res_data):
    for row in res_data.itertuples(index=False):
        tag = row.seq_id.split("|")[1]

        candidate_groups[tag].append(
            Candidate(
                id=row.candidate_id,
                start=int(row.start),
                end=int(row.end)
            )
        )


def gt_event_count(vcf):
    vcf.reset()

    for record in vcf:
        tag = record.chrom.split("|")[1]

        truth_groups[tag].append(
            Truth(
                id=len(truth_groups[tag]),
                pos=record.pos,
            )
        )
def offset_calculate(truth_groups:dict,vcf):
    vcf.reset()
    count_snp=0
    count_ins=0
    count_del=0
    count_all=0
    offset_ins=0
    offset_del=0
    offset_all=0
    for record in vcf:
        tag = record.chrom.split("|")[1]
        if tag=="SNP":
            truth_groups["SNP"][count_snp].start=record.pos
            truth_groups["SNP"][count_snp].end=record.pos
            count_snp+=1
        if tag=="INS":
            svlen=record.info.get("SVLEN")[0]
            truth_groups["INS"][count_ins].start=truth_groups["INS"][count_ins].pos+1+offset_ins
            truth_groups["INS"][count_ins].end=truth_groups["INS"][count_ins].start+svlen-1
            offset_ins+=svlen
            count_ins+=1
        elif tag=="DEL":
            end_pos=record.stop
            svlen=record.info.get("SVLEN")[0]
            truth_groups["DEL"][count_del].start=truth_groups["DEL"][count_del].pos+offset_del
            truth_groups["DEL"][count_del].end=end_pos-svlen+1+offset_del
            offset_del-=svlen
            count_del+=1
        elif tag=="ALL":
            svtype=record.info.get("SVTYPE")
            if svtype == None:
                truth_groups["ALL"][count_all].start=record.pos+offset_all
                truth_groups["ALL"][count_all].end=record.pos+offset_all
                count_all+=1
            elif svtype =="INS":
                svlen=record.info.get("SVLEN")[0]
                truth_groups["ALL"][count_all].start=truth_groups["ALL"][count_all].pos+1+offset_all
                truth_groups["ALL"][count_all].end=truth_groups["ALL"][count_all].start+svlen-1
                offset_all+=svlen
                count_all+=1
            elif svtype =="DEL":
                end_pos=record.stop
                svlen=record.info.get("SVLEN")[0]
                truth_groups["ALL"][count_all].start=truth_groups["ALL"][count_all].pos+offset_all
                truth_groups["ALL"][count_all].end=end_pos-svlen+1+offset_all
                offset_all-=svlen
                count_all+=1

def matching(truth,candidate,tolerance):
    candidate_tp=0
    truth_tp=0
    for value in candidate:
        left=0
        right=len(truth)
        start_target=value.start-tolerance
        if start_target < 0: start_target=0
        end_target=value.end+tolerance
        while left < right:
            mid = (left + right) // 2
            if truth[mid].start < start_target:
                left = mid + 1
            else:
                right = mid
        final_start=left
        left=0
        right=len(truth)
        while left < right:
            mid = (left + right) // 2
            if truth[mid].end <= end_target:
                left = mid + 1
            else:
                right = mid
        final_end=left
        if final_start < final_end:
            candidate_tp+=1
            for i in range(final_start,final_end):
                if not truth[i].matched_candidates:
                    truth_tp += 1
                truth[i].matched_candidates.append(value.id)
    result["candidate_tp"]=candidate_tp
    result["truth_tp"]=truth_tp

def main():
    parser = argparse.ArgumentParser(
        description="Read result tsv file, ground truth vcf file and metadata json file."
    )
    parser.add_argument("--res", "--result", dest="res", default="result/result.tsv",
                        help="result TSV file (default: result/result.tsv)")
    parser.add_argument("--gt", "--ground_truth", dest="gt", default="data/MS_output/mutant_ms.vcf",
                        help="ground truth VCF file (default: data/MS_output/mutant_ms.vcf)")
    parser.add_argument("--md", "--metadata", dest="md", default="result/metadata.json",
                        help="metadata json file (default: result/metadata.json)")
    args = parser.parse_args()
    res_path=project_path(args.res)
    gt_path=project_path(args.gt)
    md_path=project_path(args.md)
    res_data=pd.read_csv(res_path,sep="\t")
    with open(md_path,encoding="utf-8") as f:
        md_data=json.load(f)
    kmer_length=md_data["k-mer_length"]
    with pysam.VariantFile(gt_path) as vcf:
        first_record=next(iter(vcf),None)
        if first_record is not None and not res_data.empty :
            if first_record.chrom != res_data.iloc[0]["seq_id"]:
                raise ValueError("result tsv file and ground truth file doesn't match")
        candidate_count(res_data)
        gt_event_count(vcf)
        offset_calculate(truth_groups,vcf)
        for tag in CHROM_TAG:
            for i in TOLERANCE_RATE:
                tolerance=math.floor(i*kmer_length)
                result["candidate_total"]=len(candidate_groups[tag])
                result["truth_total"]=len(truth_groups[tag])
                matching(truth_groups[tag],candidate_groups[tag],tolerance)
                result["candidate_fp"]=result["candidate_total"]-result["candidate_tp"]
                result["truth_fn"]=result["truth_total"]-result["truth_tp"]
                print(tag)
                print("tolerance: ",tolerance)
                print("candidate_tp: ",result["candidate_tp"],"truth_tp: ",result["truth_tp"],"truth_total",result["truth_total"])
                result["candidate_fp"]=result["candidate_total"]-result["candidate_tp"]
                result["precision"]=result["candidate_tp"]/result["candidate_total"]
                result["recall"]=result["truth_tp"]/result["truth_total"]
                p=result["precision"]
                r=result["recall"]
                result["f1"]=2*p*r/(p+r) if (p+r)>0 else 0
                for data_arr in truth_groups.values():
                    for data in data_arr:
                        data.matched_candidates.clear()



if __name__ == "__main__":
    main()