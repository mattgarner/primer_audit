'''

Audit of Gemini confirmation primer failures

author@matt.j.garner@gmail.com
date@260816

Ugly code which works, but needs refactoring
'''

import csv
import sys
from datetime import date
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
#from Bio.Blast import NCBIWWW  # Enable to use BLAST_online()

def main(csv_filepath, start_date, end_date):
    
    start_date = format_date(start_date)
    print "\n\tStart date:\t", start_date
    end_date = format_date(end_date)
    print "\t  End date:\t", end_date
    
    with open(csv_filepath, 'rb') as csvfile:
        csv_reader = csv.DictReader(csvfile)
        data = extract_data(csv_reader)
        redesigns = identify_redesigned_primers(data, start_date, end_date)
        #for key in redesigns.keys():
        #    generate_primer_fasta(data, key)
        #fasta_path = evaluate_redesigns(data, redesigns)
        #fasta_path = generate_primer_fasta(data, redesigns)
        evaluate_primers(data, redesigns)


def format_date(slash_delimited_date):
    '''
    Convert a dd/mm/yyyy string into a yy/mm/dd datetime object
    '''
    try:
        date_component_order = [2,1,0]
        date_components = [int(i) for i in slash_delimited_date.split("/")]
        for number in date_components:
            assert number > 0, "Error %s !> 0" % number
        output_datetime = date( date_components[date_component_order[0]], date_components[date_component_order[1]], date_components[date_component_order[2]] )
        return output_datetime
    except:
        None


def extract_data(csv_reader):
    data_dict = {}
    try:
        for index, row in enumerate(csv_reader):
            data_dict[index+2] = row
    except csv.Error as e:
        sys.exit('file %s, line %d: %s' % (filename, csv_reader.line_num, e))    
    
    return data_dict
    
    
def identify_redesigned_primers(data_dict, start_date, end_date):
    
    records = data_dict.keys()
    total_num_primer_pairs = len(records)
    subset_num_primer_pairs = 0
    redesigned = {}
    print "\n      Primer pairs:\t", total_num_primer_pairs
    
    # for each primer record
    while records:
        record1 = records[0]
        #print "\nRecord1:", record1,
        record_date = format_date(data_dict[record1]["Date"])
        
        # filter records to range of interest
        if record_date and (record_date <= end_date) and (record_date >= start_date):
            subset_num_primer_pairs += 1
            
            # search for a matching record in remaining records
            for index, record2 in enumerate(records[1:]):
                # records with same gpos and designed after record1
                if (data_dict[record2]["Co-ordinate"] == data_dict[record1]["Co-ordinate"]) and (record_date < format_date(data_dict[record2]["Date"])):
                
                    # record with different primer(s) suggesting failure of 1st pair
                    # first remove tags from primer seqs
                    #print "\n"
                    r1_primer_F = trim_primer_tags(data_dict[record1]["Fseq"])
                    r1_primer_R = trim_primer_tags(data_dict[record1]["Rseq"])
                    r2_primer_F = trim_primer_tags(data_dict[record2]["Fseq"])
                    r2_primer_R = trim_primer_tags(data_dict[record2]["Rseq"])
                    #print "1R: ", r1_primer_R
                    #print "2R: ", r2_primer_R
                    
                    # then check for diff F or R seq, indicating failure and redesign
                    if (r2_primer_F != r1_primer_F) or (r2_primer_R != r1_primer_R):
                        # Save record pairings for redesigns
                        try:
                            redesigned[record1].append(record2)
                        except KeyError:
                            redesigned[record1] = [record2]
                             
                        # remove the match from records so we don't match the pair twice
                        #print "\ndeleting r2", records[index+1]
                        del records[index+1]
        #print "deleting ", records[0]
        del records[0]
        
    print " Within date range:\t", subset_num_primer_pairs
    count_redesigns = sum(len(redesigned[key]) for key in redesigned.keys())  
    print "Primers redesigned:\t", len(redesigned.keys())
    print "  Failed redesigns:\t", count_redesigns - len(redesigned.keys())
    print "   Total redesigns:\t", count_redesigns
    return redesigned


#===============================================================================
# def generate_primer_fasta(data, subset):
# 
#     fasta_filepath = "./fasta.fa"
#     for row in subset.keys():
#         #row = str(key+2)
#         write_fasta(row + "_F", data[key]["Fseq"], fasta_filepath)
#         write_fasta(row + "_R", data[key]["Rseq"], fasta_filepath)
#                         
#         for row in subset[key][:-1]:    
#             #row = value + 2
#             write_fasta(str(row) + "_F", data[value]["Fseq"], fasta_filepath)
#             write_fasta(str(row) + "_R", data[value]["Rseq"], fasta_filepath)
#         
#         #row = str(subset[key][-1]+2)
#         write_fasta(row + "_F", data[subset[key][-1]]["Fseq"], fasta_filepath)
#         write_fasta(row + "_R", data[subset[key][-1]]["Rseq"], fasta_filepath)
#         
#     return fasta_filepath
#===============================================================================


def generate_primer_fasta(data, records, fasta_filepath = "./fasta.fa"):

    primers = {"forward":{"suffix":"_F", 
                          "tag":"Fseq"}, 
               "reverse":{"suffix":"_R", 
                          "tag":"Rseq"}
               }

    with open(fasta_filepath, "w") as fasta_fh:    
        for record in records:
            for primer in primers.keys():
                header = ">" + str(record) + primers[primer]["suffix"]+"\n"
                primer_seq = trim_primer_tags(data[record][primers[primer]["tag"]]) + "\n"
                
                fasta_fh.write(header)
                fasta_fh.write(primer_seq)

        #write_fasta(str(record) + "_F", data[record]["Fseq"], fasta_filepath)
        #write_fasta(str(record) + "_R", data[record]["Rseq"], fasta_filepath)

    return fasta_filepath


def evaluate_redesigns(data, redesigns):
    #print redesigns
    print "\n\n###Redesigns###"
    fasta_filepath = "./fasta.fa"
    i = 1
    for key in redesigns.keys():
        print "\n########################################################\n\n",i
        print "Original design: "
        row = str(key+2)
        print row, data[key]["Date"]
        print "F: ", data[key]["Fseq"]
        write_fasta(row + "_F", data[key]["Fseq"], fasta_filepath)
        print "R: ", data[key]["Rseq"]
        write_fasta(row + "_R", data[key]["Rseq"], fasta_filepath)
                        
        print "\nFailed redesigns: "
        for value in redesigns[key][:-1]:    
            row = value + 2
            print row, data[value]["Date"]
            print "F: ", data[value]["Fseq"]
            write_fasta(str(row) + "_F", data[value]["Fseq"], fasta_filepath)
            print "R: ", data[value]["Rseq"]
            write_fasta(str(row) + "_R", data[value]["Rseq"], fasta_filepath)
        
        print "\nSuccessful redesign:"
        row = str(redesigns[key][-1]+2)
        print row, data[redesigns[key][-1]]["Date"]
        print "F: ", data[redesigns[key][-1]]["Fseq"]
        write_fasta(row + "_F", data[redesigns[key][-1]]["Fseq"], fasta_filepath)
        print "R: ", data[redesigns[key][-1]]["Rseq"]
        write_fasta(row + "_R", data[redesigns[key][-1]]["Rseq"], fasta_filepath)
        
        extra_time = format_date(data[redesigns[key][-1]]["Date"]) - format_date(data[key]["Date"])
        print "\nAdditional time taken: %s days" % extra_time.days
    
        i += 1
    return fasta_filepath


def evaluate_primers(data, redesigns):
    
    for key in redesigns.keys():
        print "\n\n\n#####################################\n\n\n"
        print key, "FAILED FIRST DESIGN"
        evaluate_primer_pair(data, key)
        for value in redesigns[key][:-1]:
            print "\n\n"
            print "%s FAILED REDESIGN OF %s" % (value, key)
            evaluate_primer_pair(data, value)
        print "\n\n"
        print "%s SUCCESSFUL REDESIGN OF %s" % (redesigns[key][-1], key)
        evaluate_primer_pair(data, redesigns[key][-1])
        
        delay = format_date(data[redesigns[key][-1]]["Date"]) - format_date(data[key]["Date"])
        print "\nAdditional time taken: %s days" % delay.days
        

def evaluate_primer_pair(data, record):
        
        print "Date:", data[record]["Date"]
        print "Target:", data[record]["Co-ordinate"]
        fasta_filepath = generate_primer_fasta(data, [record])

        # Blast primers
        blast_xml_results_path = BLAST_local(fasta_filepath)
        blast_xml_results_path = "blast_results.xml"
        blast_results = open(blast_xml_results_path)

        # Parse output
        records = NCBIXML.parse(blast_results)

        # For each query sequence (i.e. primer)
        forward_blast = next(records)
        reverse_blast = next(records)
        
        best_forward_alignment = None
        best_reverse_alignment = None
        
        try:
            best_forward_alignment = find_best_alignments(forward_blast)[0][0]
        except:
            print "No best alignment for F primer"
        try:
            best_reverse_alignment = find_best_alignments(reverse_blast)[0][0]
        except:
            print "No best alignment for R primer"
        
        print "\n\nBEST MAPPING EVALUATION"
        if best_forward_alignment and best_reverse_alignment:
            evaluate_position(best_forward_alignment, best_reverse_alignment, data[record]["Co-ordinate"])
        

def get_strand(hsp):
    if (hsp.sbjct_start < hsp.sbjct_end):
        strand = "+"
    else:
        strand = "-"
    return strand


def evaluate_position(forward_alignment, reverse_alignment, amplicon_target):
    # Possible metrics:
        # are the primers in the correct region
        # what is the target
        # do primers map within x bp
        # distance between coords
        # target lies between coords 
        # opposite target strands
        # primer dimers
            # alignment score
        # specificity - edit distance to next best match
        # Tm diff between each in pair
        # SNPs - should be checked during design

    
    target_chr, target_coord = amplicon_target.split(":")
    print "Intended target:", amplicon_target
    
    forward_mapping, forward_hsp = forward_alignment
    reverse_mapping, reverse_hsp = reverse_alignment    
    
    # Check chr mapping is correct
    forward_mapping_chr = forward_mapping.split(" ")[0]
    reverse_mapping_chr = reverse_mapping.split(" ")[0]
    
    required_strands = ["+","-"]
    forward_primer_strand = get_strand(forward_hsp)
    reverse_primer_strand = get_strand(reverse_hsp)
        
    for strand in required_strands:
        if not strand in [forward_primer_strand, reverse_primer_strand]:
            print "No primer on %s strand!" % strand
            break
    else:
        print "Primers on opposite strands"
    
    if forward_mapping_chr == reverse_mapping_chr:
        same_mapping_chr = True
        print "Both primers map to same chr: %s" % forward_mapping_chr
        if forward_mapping_chr == target_chr:
            print "Primers map to the target chr: %s" % forward_mapping_chr
        else:
            print "Primers map to the WRONG chr: %s" % forward_mapping_chr
    
    else:
        same_mapping_chr = False
        print "Primers map to different chrs! F: %s R: %s" % (forward_mapping_chr, reverse_mapping_chr) 
    
    # Check target coord is flanked by primers of opposite strands
    amplicon_start = min(forward_hsp.sbjct_start, forward_hsp.sbjct_end, reverse_hsp.sbjct_start, reverse_hsp.sbjct_end)
    amplicon_end = max(forward_hsp.sbjct_start, forward_hsp.sbjct_end, reverse_hsp.sbjct_start, reverse_hsp.sbjct_end)
    amplicon_length = amplicon_end-amplicon_start+1
    
    lower_target_coord = int(target_coord.split("-")[0])
    upper_target_coord = int(target_coord.split("-")[-1])
    start_to_target_distance = lower_target_coord - amplicon_start
    target_to_end = amplicon_end - upper_target_coord
    target_length = upper_target_coord - lower_target_coord + 1
    
    if (amplicon_start < lower_target_coord) and (amplicon_end > upper_target_coord):
        target_within_amplicon = True
    else:
        target_within_amplicon = False
    
    print "F primer mapping coords:", forward_hsp.sbjct_start, forward_hsp.sbjct_end
    print "R primer mapping coords:", reverse_hsp.sbjct_start, reverse_hsp.sbjct_end
    print "Amplicon range: %s - %s" % (amplicon_start, amplicon_end)
    print "Amplicon length:", amplicon_length
    print "Target within amplicon:", target_within_amplicon
    
        
    flanking_target = True
    amplicon_length = 0
    in_range = True
    opp_strands = True
    return True

        
def find_best_alignments(blast_record):
        
    best_matches = None
    second_best_matches = None
    
    print "\n\n\n" + blast_record.query 
    print "query length", blast_record.query_length
    # For each sequence aligned to (i.e.chr)
    for alignment in blast_record.alignments:
                    
        for hsp in alignment.hsps:
            
            # No best exists yet
            if best_matches == None:
                best_matches = [[alignment.title, hsp]]
                        
            # New best score
            elif hsp.score > best_matches[-1][1].score:
                second_best_matches = [best_matches]
                best_matches = [[alignment.title, hsp]]
            
            # Equal best score
            elif hsp.score == best_matches[-1][1].score:
                best_matches.append([alignment.title, hsp])
            
            # No second best exists yet
            elif second_best_matches == None:
                second_best_matches = [[alignment.title, hsp]]
            
            # Equal second score
            elif hsp.score == second_best_matches[-1][1].score:
                second_best_matches.append([alignment.title, hsp])
    
    print "Best:"
    if best_matches:
        for item in best_matches:
            hsp = item[1]
            print item[0], "\t", int(hsp.score), "\\", blast_record.query_length
    
            print('****Alignment****')
            print('aligned length:', hsp.align_length)
            print('qstart:', hsp.query_start)
            print('qend:', hsp.query_end)
            print('sstart:', hsp.sbjct_start)
            print('send:', hsp.sbjct_end)
            print('score:', hsp.score)
            print('gaps:', hsp.gaps)
            print('e value:', hsp.expect)
            if (hsp.sbjct_start < hsp.sbjct_end):# and (query.query[-1] == "F")) or\
               #((hsp.sbjct_start > hsp.sbjct_end) and (query.query[-1] == "R")):
                strand = "+"
            else:
                strand = "-"
            print('strand:', strand)
            print('frame:', hsp.frame)
            print(hsp.query)
            print(hsp.match)
            print(hsp.sbjct)
    
        
        
        print "Next best matches:"
        for item in second_best_matches:
            record = item[0]
            hsp = item[1]
            print record, "\t", int(hsp.score), "\\", blast_record.query_length
        
        if len(best_matches) > 1:
            print "Multiple best matches detected!"
        else:
            print "Edit distance:", int(best_matches[-1][1].score - second_best_matches[-1][1].score)    
        
        return [best_matches, second_best_matches]
    return None


def BLAST_online(forward, reverse):
    # Not used
    print "\nPerforming BLAST search for %s" % sequence 
    result_handle = NCBIWWW.qblast("blastn-short", "refseq_genomic", sequence, entrez_query="txid9606[ORGN]", hitlist_size=5)
    save_file = open("my_blast.xml", "w")
    save_file.write(result_handle.read())
    save_file.close()
    result_handle.close()
    print "search complete!"


def BLAST_local(fasta, output_filepath = "./blast_results.xml"): 
    blastn_cline = NcbiblastnCommandline(query="fasta.fa", task='blastn-short', db="~/primer_audit/human_1KG_v37_BLASTdb/human_1KG_v37_BLASTdb", out=output_filepath, outfmt=5)# )
    stdout, stderr = blastn_cline()
    assert stdout == ''
    assert stderr == ''
    return output_filepath


def write_fasta(seq_id, seq, fasta_filepath):
    header = ">" + str(seq_id)
    with open(fasta_filepath, "a") as fasta_fh:
        fasta_fh.write(header+"\n")
        fasta_fh.write(trim_primer_tags(seq)+"\n")


def trim_primer_tags(seq):
    return "".join([char for char in seq if char.isupper()])


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])