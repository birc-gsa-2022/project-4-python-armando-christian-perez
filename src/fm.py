import argparse
import sys

def main():
    argparser = argparse.ArgumentParser(
        description="FM-index exact pattern matching",
        usage="\n\tfm -p genome\n\tfm genome reads"
    )
    argparser.add_argument(
        "-p", action="store_true",
        help="preprocess the genome."
    )
    argparser.add_argument(
        "genome",
        help="Simple-FASTA file containing the genome.",
        type=argparse.FileType('r')
    )
    argparser.add_argument(
        "reads", nargs="?",
        help="Simple-FASTQ file containing the reads.",
        type=argparse.FileType('r')
    )
    args = argparser.parse_args()
    name = args.genome.name.replace(".fa", "")
    if args.p:
        create_BWT_structures(args.genome, name)
        
    else:
        fasta_dict, fastq_dict = fasta_translator(args.genome), fastq_translator(args.reads)
        SAM = matches_to_SAM(fasta_dict, fastq_dict, name)
        print_SAM(SAM)
        # here we need the optional argument reads
        if args.reads is None:
            argparser.print_help()
            sys.exit(1)

def fasta_translator(input_file):
    output_dict = {}
    start = True
    for i in input_file:
        i = i.strip()
        if len(i) == 0:
            continue
        if i[0] == ">":
            if  not start:
                output_dict[name] = seq + "$"
            name = i[1:].strip()
            seq = ""
            if start:
                start = False
        elif start:
            continue
        else:
            seq += i
    output_dict[name] = seq + "$"
    return output_dict

def fastq_translator(input_file):
    output_dict = {}
    for i in input_file:
        if i[0] == "@":
            name = i[1:].strip()
        else:
            seq = i.strip()
            output_dict[name] = seq
    return(output_dict)

def compute_buckets(SA, reference):
    bucket_len = {}
    for i in SA:
        if reference[i - 1] in bucket_len:
            bucket_len[reference[i - 1]] += 1
        else:
            bucket_len[reference[i - 1]] = 1
    bucket_len = dict(sorted(bucket_len.items()))
    buckets = {}
    len_traversed = 0

    for key in bucket_len:
        temp = len_traversed
        len_traversed += bucket_len[key]
        buckets[key] = temp

    return buckets

def create_tracker(buckets, SA):
    tracker = {}
    for i in buckets:
        tracker[i] = [None for i in range(len(SA) + 1)]
    return tracker

def create_ranks(buckets, SA, reference):
    tracker = create_tracker(buckets, SA)
    for char in tracker:
        tracker[char][0] = 0
    for i in range(1, len(SA) + 1):
        for char in tracker:
            tracker[char][i] = tracker[char][i - 1]
        prev_bwt_char = reference[SA[i - 1] - 1]
        tracker[prev_bwt_char][i] += 1
    return tracker

def SA_naive(string):
    SA = [i for i in range(len(string))]
    SA = [string[SA[int(i)]:] for i in SA]
    SA.sort()
    SA = [len(string) - len(i) for i in SA]
    return(SA)

def create_BWT_structures(ref_file, name):
    fasta_dict = fasta_translator(ref_file)

    SA_file = open(name + r"_SA.txt", "w")
    bucket_file = open(name + r"_buckets.txt", "w")
    rank_file = open(name + r"_rank.txt", "w")

    for key in fasta_dict:
        SA = SA_naive(fasta_dict[key])
        SA_print = ",".join(str(i) for i in SA)
        SA_file.write(SA_print)
        buckets = compute_buckets(SA, fasta_dict[key])
        for bucket in buckets:
            bucket_file.write(bucket + " " + str(buckets[bucket]) + "\t")
        rank = create_ranks(buckets, SA, fasta_dict[key])
        rank_chars = list(rank.keys())
        rank_file.write(",".join(rank_chars))
        for i in range(len(SA) + 1):
            rank_file.write("\t")
            row = [str(rank[char][i]) for char in rank_chars]
            rank_file.write(",".join(row))
        rank_file.write("\n")
        SA_file.write("\n")
        bucket_file.write("\n")

    SA_file.close()
    bucket_file.close()
    rank_file.close()

def read_SA(line):
    SA = line.strip().split(",")
    SA = [int(sa) for sa in SA]
    return SA

def read_bucket(line):
    buckets = {}
    ite = line.strip().split("\t")
    for i in ite:
        info = i.split(" ")
        buckets[info[0]] = int(info[1])
    return buckets

def read_rank(line):
    ite = line.strip().split("\t")
    keys = ite[0].split(",")
    ranks = {key: [None for i in range(len(ite) - 1)] for key in keys}
    for i in range(1, len(ite)):
        row = [int(c) for c in ite[i].split(",")]
        for j, key in enumerate(keys):
            ranks[key][i - 1] = row[j]
    return ranks

def matches_to_SAM(fasta_dict, fastq_dict, name):
    read_name = []
    reference_name = []
    match_index = []
    CIGARS = []
    match_string = []

    keys = [i for i in fasta_dict]

    SA_file = open(name + r"_SA.txt", "r")
    bucket_file = open(name + r"_buckets.txt", "r")
    rank_file = open(name + r"_rank.txt", "r")
    for i, (SA_line, bucket_line, rank_line) in enumerate(zip(SA_file, bucket_file, rank_file)):
        SA = read_SA(SA_line)
        buckets = read_bucket(bucket_line)
        ranks = read_rank(rank_line)
        for read_key in fastq_dict:
            matches = FM_match(fasta_dict[keys[i]], fastq_dict[read_key], SA, buckets, ranks)
            for match in matches:
                reference_name.append(keys[i])
                read_name.append(read_key)
                match_index.append(match + 1)
                CIGARS.append(str(len(fastq_dict[read_key])) + "M")
                match_string.append(fastq_dict[read_key])
    SA_file.close()
    bucket_file.close()
    rank_file.close()
    output = (read_name, reference_name, match_index, CIGARS, match_string)
    return output

def FM_match(reference, read, SA, buckets, ranks):
    upper = 0
    lower = len(SA)
    for char in reversed(read):
        if upper == lower:
            return []
        upper = buckets[char] + ranks[char][upper] # Will break if read char is not in genome
        lower = buckets[char] + ranks[char][lower]
    matches = SA[upper:lower]
    return matches

def print_SAM(SAM):
    for i in range(len(SAM[0])):
        sys.stdout.write(SAM[0][i] + "\t" + SAM[1][i] + "\t" + str(SAM[2][i]) + "\t" + SAM[3][i] + "\t" + SAM[4][i] + "\n")

if __name__ == '__main__':
    main()