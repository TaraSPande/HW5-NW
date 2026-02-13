# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/Tursiops_truncatus_BRD2.fa")

    # Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix

    # Print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix

    nw = NeedlemanWunsch(sub_matrix_file="./substitution_matrices/BLOSUM50.mat", gap_open=-10, gap_extend=-1)

    alignment = {}
    for header, seq in {gg_header: gg_seq, mm_header: mm_seq, br_header: br_seq, tt_header: tt_seq}.items(): #for each species
        score, seqA_align, seqB_align = nw.align(hs_seq, seq)   #align human to current species
        alignment[header] = score   #get score

    sorted_species = sorted(alignment.items(), key=lambda x: x[1], reverse=True)    #sort from best to worst

    for header, score in sorted_species:    #print out species in order with score
        species = header.split("OS=")[1].split(" ")[0]
        print(f"Human—{species}: {score}")

if __name__ == "__main__":
    main()
