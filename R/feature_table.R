#' Build table of binary features
#'
#' @param con configuration
#' @param pair a list of alignment result
#' @export
#' @return a data frame of feature table
buildFeatureTable <- function(con, pair){

  genes = pair$gene_tag
  gene2genome = pair$gene2genome.pval
  gaps = pair$gaps
  indels = pair$indels
  aa2aa = pair$aa2aa.pval
  aa2orf = pair$aa2orf.pval
  aa2transorf = pair$aa2transorf.pval
  synder_summary = pair$summary_synder_flags
  synder_out = pair$synder_out
  skipped = pair$skipped
  genome_hits = pair$gene2genome

  "Set p-value cutoffs for each set of alignments (aa-vs-aa, aa-vs-orf,
    aa-vs-mRNA, gene-vs-genome). The p-value threshold set in the configuration
    is the overall desired p-value (defaulting to 0.05)."

  p2p_cutoff <- con@alignment@pvalue@prot2prot
  p2a_cutoff <- con@alignment@pvalue@prot2allorf
  d2d_cutoff <- con@alignment@pvalue@dna2dna
  p2t_cutoff <- con@alignment@pvalue@prot2transorf

  p2p_alncut <- con@alignment@alnrate@prot2prot
  p2a_alncut <- con@alignment@alnrate@prot2allorf
  d2d_alncut <- con@alignment@alnrate@dna2dna
  p2t_alncut <- con@alignment@alnrate@prot2transorf


  met <- merge(
    gene2genome,
    genome_hits$map[, c("query", "cds_match", "exon_match", "mrna_match")],
    by="query"
  ) %>%
    dplyr::select(-target, -pval) %>%
    dplyr::filter(pval.adj < d2d_cutoff & alnrate > d2d_alncut)

  feature_table <- data.frame(
    seqid = genes,
    # matches somewhere in at least one search interval
    nuc = genes %in% met$query,
    # matches CDS in at least one search interval
    cds = genes %in% subset(met, cds_match)$query,
    # matches transcript in at least one search interval
    rna = genes %in% subset(met, mrna_match)$query,
    # matches exon in at least one search interval
    exo = genes %in% subset(met, exon_match)$query,
    # at least search interval overlaps a N-string
    nst = genes %in% gaps$query,
    # number of confirmed indels (based on search interval size)
    ind = genes %in% indels,
    # the query has an ortholog in the target
    gen = genes %in% subset(aa2aa, pval.adj < p2p_cutoff & alnrate > p2p_alncut)$query,
    # ORF match in SI
    orf = genes %in% subset(aa2orf, pval.adj < p2a_cutoff & alnrate > p2a_alncut)$query,
    # ORF match to spliced transcript (possibly multi-exonic)
    trn = genes %in% subset(aa2transorf, pval.adj < p2t_cutoff & alnrate > p2t_alncut)$query,
    # synteny is scrambled
    scr = genes %in% subset(synder_summary, incoherent)$attr,
    # at least one search interval maps off scaffold
    una = genes %in% subset(synder_summary, unassembled)$attr,
    # search interval was not processed for technical reasons (e.g. too big)
    tec = genes %in% skipped,
    stringsAsFactors = FALSE
  )
  feature_table
}


#' Build label tree for homolog decision
#'
#' @param feats a data frame of feature table
#' @param con configuration
#' @export
#' @return a label tree for feature table
buildLabelsTree <- function(feats, con){

  "Merge all features into the feature decision tree. Each node contains a
  logical vector specifying which queries are members of the node. The names
  in the decision tree must match the names in the feature table."

  root <- con@decision_tree %>%
    data.tree::as.Node(replaceUnderscores=FALSE)
  classify <- function(node, membership=logical(0)){
    if(length(membership) == 0){
      membership <- rep(TRUE, nrow(feats))
    }
    node$membership <- membership
    node$N <- sum(membership)
    if(node$name %in% names(feats)){
      kids <- node$children
      if(length(kids) == 2){
        yes <-  feats[[node$name]] & membership
        no  <- !feats[[node$name]] & membership
        classify(kids[[1]], yes)
        classify(kids[[2]], no)
      } else if (length(kids) != 0) {
        stop("Nodes have either 0 or 2 children")
      }
    } else if(node$isRoot){
      classify(node$children[[1]], membership)
    }
  }
  classify(root)
  root
}


#' Match label to feature table
#'
#' @param root a label tree for feature table
#' @param feats a data frame of feature table
#' @export
#' @return a data frame combine homolog class label and feature table
labelTreeToTable <- function(root, feats){

  "
  Build a dataframe for each node
  For example:
          seqid primary secondary
  1 AT1G29418.1   ORFic        O1
  2 AT3G15909.1   ORFic        O1
  Then bind the rows of the dataframes for each node.
  "

  toTable <- function(node) {
    if(!is.null(node$N) && node$N > 0){
      seqid     <- feats$seqid[node$membership]
      primary   <- node$primary
      secondary <- node$secondary
      # If there are no members in the class
    } else {
      seqid     <- character(0)
      primary   <- character(0)
      secondary <- character(0)
    }
    data.frame(
      seqid     = seqid,
      primary   = primary,
      secondary = secondary,
      stringsAsFactors=FALSE
    )
  }
  # Get a table for each class
  # --- NOTE: Without simplify=FALSE something truly dreadful happens.
  # The proper output (where a,b,c ... are dataframe columns):
  #   <Node_1> a b c
  #   <Node_2> d e f
  #   <Node_3> g h i
  #   <Node_4> j k l
  # The automagically borken output
  #   <Node_1> a
  #   <Node_2> b
  #   <Node_3> c
  #   <Node_4> d
  # The columns are recast as vectors, fed to the wrong children, and the
  # leftovers are tossed.
  root$Get(toTable, filterFun = data.tree::isLeaf, simplify=FALSE) %>%
    # Bind all tables into one
    do.call(what=rbind) %>%
    # Remove rownames
    set_rownames(NULL)
}


#' Summarize homolog class for the list of genes from focal species
#'
#' @param con configuration
#' @param pair a list of alignment result
#' @param tarname a string of target species name
#' @export
#' @return a data frame of feature table
merge_feature_table <- function(con, pair, tarname){

  name_conversion <-  c(O1="A_gen", O2="A_trn", O3="A_orf",
                        N1="N_cds", N2="N_exo", N3="N_rna", N4="N_dna",
                        U2="U_ind", U5="U_scr", U6="U_unk", U1="U_una", U3="U_nst", U7="U_tec")

  feature_table <- buildFeatureTable(con, pair)
  root <- buildLabelsTree(feature_table,con)
  labels <- labelTreeToTable(root,feature_table) %>%
    dplyr::select(seqid, homology_class=secondary) %>% dplyr::mutate(homology_class = name_conversion[homology_class])
  final <- merge(feature_table, labels, by='seqid')
  #final$homology_class[final$homology_class=="A_trn" & final$cds==T] <- "N_cds"
  final$target_species <- tarname
  final
}


