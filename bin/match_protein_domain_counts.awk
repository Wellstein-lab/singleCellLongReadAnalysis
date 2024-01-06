# program:  match_protein_domain_counts.awk
# purpose:  To create a matrix of the experiments by the protein and the domains
#           and the read counts
#
# input:  the output from the countDomainPerSequence.sh
#
# output: matrix where the counts of all the experiment read counts  are the columns
#         rows are the defining series of protein, domain, aa-position and aa-sequence)
#
# Steps:  #
# Read file B, each of the new sorted arrays and store the matching columns in an associative array
# an associative array is a key-value pair
# where the key is the unique identifier that yields the value
# in our case the values are the read counts  which are unique to this
# sample
#
# Each individual array is as follows:
#
#          col 1 - Protein
#          col 2 - Domain_Name
#          col 3 - AA_position
#          col 4 - Domain_sequence
#          col 5 - Read_counts

# A is the union file and is as follows
#
#          col 1 - Protein
#          col 2 - Domain_Name
#          col 3 - AA_position
#          col 4 - Domain_sequence
#
#
# The match will print out the content of A, the union file, and the Counts from sample B, the individual results for the experiment with read counts
#
# this ensures all the individual samples files have the same structure and IDs for making a single matrix
#
NR == FNR {
  # First file: Read and store data in an array
  key = $1 OFS $2 OFS $3 OFS $4
  value = $5
  matchArray[key] = value
  next
}

{
  # Second file: Merge data based on matching keys
  key = $1 OFS $2 OFS $3 OFS $4
  value = $5
  if (key in matchArray) {
    print key, matchArray[key], value
  }
}
