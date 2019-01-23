source("./function.coevo.R")
require( stringr )

# examine the positive site detect by MEME
anc_139 = parse_hyphy( folerdif = "./selection/out/pH5_c2344anc_139/" )
anc_123 = parse_hyphy( folerdif = "./selection/out/pH5_c2344anc_124/" )
anc_88  = parse_hyphy( folerdif = "./selection/out/pH5_c2344anc_88/" )

c2344_clade = parse_hyphy( folerdif = "./selection/out/pH5_c2344_clades/", n_sample = 4)


meme_resi_anc   <- intersect( anc_123[[1]][[4]],  anc_88[[1]][[4]] )
meme_resi_clade <- c( c2344_clade[[1]][[4]], c2344_clade[[2]][[4]], c2344_clade[[3]][[4]], c2344_clade[[4]][[4]] )


ha_num( ref_fas = "./aa_numbering/align_aa_ref.fasta", ref_csv = "./aa_numbering/ref_numbering.csv", data_pos = meme_resi_clade )
