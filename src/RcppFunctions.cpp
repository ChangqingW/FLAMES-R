#include <Rcpp.h>

#include "main-functions/flexiplex.h"
//' Rcpp port of flexiplex
//'
//' @description demultiplex reads with flexiplex, for detailed description, see
//' documentation for the original flexiplex: https://davidsongroup.github.io/flexiplex
//'
//' @param r_segments List defining the barcode structure
//' @param r_barcode_groups List defining barcode groups
//' @param max_flank_editdistance int, maximum edit distance for matching flanking sequences
//' @param reads_in Input FASTQ or FASTA file
//' @param reads_out output file for demultiplexed reads
//' @param stats_out output file for demultiplexed stats
//' @param bc_out WIP
//' @param reverseCompliment bool, whether to reverse complement the reads after demultiplexing
//' @param n_threads number of threads to be used during demultiplexing
//' @return integer return value. 0 represents normal return.
//' @export
// [[Rcpp::export]]
Rcpp::IntegerVector flexiplex(
  Rcpp::List r_segments,
  Rcpp::List r_barcode_groups,
  int max_flank_editdistance,
  Rcpp::StringVector reads_in,
  Rcpp::String reads_out,
  Rcpp::String stats_out,
  Rcpp::String bc_out, bool reverseCompliment, int n_threads) {

  return flexiplex_cpp(
    r_segments,
    r_barcode_groups,
    max_flank_editdistance,
    reads_in,
    reads_out,
    stats_out,
    bc_out, reverseCompliment, n_threads
  );
}

#include "main-functions/pileup_readid.h"
// [[Rcpp::export]]
Rcpp::NumericMatrix
variant_count_matrix(Rcpp::String bam_path, Rcpp::String seqname, int pos,
                     bool indel, bool verbose) {
  return variant_count_matrix_cpp(bam_path, seqname, pos, indel,
                                  verbose);
}
