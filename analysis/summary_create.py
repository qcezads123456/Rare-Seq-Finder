from __future__ import annotations
import argparse
import pandas as pd
from pathlib import Path
import pysam
import json
import math
from dataclasses import dataclass,field
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
            offset_del-=abs(svlen)
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
                offset_all-=abs(svlen)
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
    return candidate_tp,truth_tp

def header_write_in(summary_set:list,multimatch_set:list):
    pd.DataFrame(columns=summary_set[0]).to_csv(
        summary_set[1],
        sep="\t",
        index=False
    )
    pd.DataFrame(columns=multimatch_set[0]).to_csv(
        multimatch_set[1],
        sep="\t",
        index=False
    )

def data_write_in(data_set:list):
    pd.DataFrame([data_set[2]]).to_csv(
        data_set[1],
        sep="\t",
        header= False,  
        index=False,
        mode="a"
        )

def summary_encapsulate(tag,tolerance,summary_value:list):
    summary_value.append(len(candidate_groups[tag]))
    candidate_tp,truth_tp=matching(truth_groups[tag],candidate_groups[tag],tolerance)
    summary_value.append(candidate_tp)
    summary_value.append(summary_value[2]-summary_value[3])
    summary_value.append(len(truth_groups[tag]))
    summary_value.append(truth_tp)
    summary_value.append(summary_value[5]-summary_value[6])
    #print(tag)
    #print("tolerance: ",tolerance)
    #print("candidate_tp: ",result["candidate_tp"],"truth_tp: ",result["truth_tp"],"truth_total",result["truth_total"])
    summary_value.append(round(summary_value[3]/summary_value[2],4) if summary_value[2] > 0 else 0.0)
    summary_value.append(round(summary_value[6]/summary_value[5],4) if summary_value[5] > 0 else 0.0)
    p=summary_value[8]
    r=summary_value[9]
    summary_value.append(round(2*p*r/(p+r),4) if (p+r)>0 else 0.0)

def multimatch_encapsulate(tag,tolerance,multimatch_value:list,record):
    multimatch_value.append(record.id+1)
    multimatch_value.append(tag)
    multimatch_value.append(tolerance)
    multimatch_value.append(record.matched_candidates)

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
    with open(md_path,encoding="utf-8") as j_f:
        md_data=json.load(j_f)
    kmer_length=md_data["k-mer_length"]
    with pysam.VariantFile(gt_path) as vcf:
        first_record=next(iter(vcf),None)
        if first_record is not None and not res_data.empty :
            if first_record.chrom != res_data.iloc[0]["seq_id"]:
                raise ValueError("result tsv file and ground truth file doesn't match")
        candidate_count(res_data)
        gt_event_count(vcf)
        offset_calculate(truth_groups,vcf)
        serial_num=first_record.chrom.split("|")[0]
        summary_path=project_path(f"analysis/summary/{serial_num}_validation.tsv")
        multimatch_path=project_path(f"analysis/summary/{serial_num}_multimatch.tsv")
        summary_header=[
            "profile","tolerance length","candidate total","candidate tp",
            "candidate fp","truth total","truth tp","truth fn",
            "precision","recall","F1-score"
        ]
        multimatch_header=["event ID","profile","tolerance length","matched candidate"]
        summary_value=[]
        multimatch_value=[]
        summary_set=[summary_header,summary_path,summary_value]
        multimatch_set=[multimatch_header,multimatch_path,multimatch_value]
        header_write_in(summary_set,multimatch_set)
        for tag in CHROM_TAG:
            for i in TOLERANCE_RATE:
                tolerance=math.floor(i*kmer_length)
                summary_value.append(tag)
                summary_value.append(tolerance)
                summary_encapsulate(tag,tolerance,summary_value)
                data_write_in(summary_set)
                for record in truth_groups[tag]:
                    multimatch_encapsulate(tag,tolerance,multimatch_value,record)
                    data_write_in(multimatch_set)
                    multimatch_value.clear()
                for data_arr in truth_groups.values():
                    for data in data_arr:
                        data.matched_candidates.clear()
                summary_value.clear()
                multimatch_value.clear()


if __name__ == "__main__":
    main()