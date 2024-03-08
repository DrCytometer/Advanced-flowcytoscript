# Liston/Dooley lab flowcytoscript

# Non-interactive version for users comfortable with R.
# More control over output and analysis parameters is provided
# via flowcytoscript_parameter file.


## Start-up ==================
library( digest )
require( dunn.test )
library( ggplot2 )
library( ggridges )
library( ggrepel )
library( ggtext )
library( RColorBrewer )
library( Rtsne )
library( uwot )
library( dplyr )
library( tidyr )
library( FastPG )
require( RcppHNSW )
library( parallel )
library( readxl )
library( emmeans )
library( ConsensusClusterPlus )
library( EmbedSOM )
library( flowCore )
library( flowWorkspace )
library( data.table )
library( scattermore )


# source parameters and functions---------------

# source parameter file

param.filename <- "./flowcytoscript_parameter.r"
source( param.filename )

source( file.path( src.dir, "ce_diff_test.r" ) )
source( file.path( src.dir, "ce_diff_test_tsne.r" ) )
source( file.path( src.dir, "ce_diff_test_umap.r" ) )
source( file.path( src.dir, "plot_cluster_dmrd.r" ) )
source( file.path( src.dir, "plot_all_dmrd_figures.r" ) )
source( file.path( src.dir, "plot_dimensionality_reduction.r" ) )
source( file.path( src.dir, "ggplot_cluster_dmrd.r" ) )
source( file.path( src.dir, "ggplot_all_dmrd_figures.r" ) )
source( file.path( src.dir, "ggplot_dimensionality_reduction.r" ) )
source( file.path( src.dir, "plot_changed_regions.r" ) )


# stopifnot checks

stopifnot( names( fcs.channel.label ) == fcs.channel )

stopifnot( names( fcs.condition.label ) == fcs.condition )
stopifnot( names( fcs.condition.color ) == fcs.condition )
stopifnot( names( fcs.condition.line.type ) == fcs.condition )

stopifnot( names( fcs.cluster.label ) == fcs.cluster )
stopifnot( names( fcs.cluster.color ) == fcs.cluster )
stopifnot( names( fcs.cluster.line.type ) == fcs.cluster )

stopifnot( sum( is.null( fcs.dmrd.data.sample.n ), 
                is.null( fcs.dmrd.data.sample.n.per.condition ), 
                is.null( fcs.dmrd.data.sample.n.per.sample ) ) >= 2 )

# create dirs

figure.dir <- c( 
  fcs.ce.diff.tsne.figure.dir, 
  fcs.ce.diff.umap.figure.dir, 
  fcs.density.figure.dir, 
  fcs.heatmap.figure.dir, 
  fcs.histogram.figure.dir, 
  fcs.tsne.figure.dir, 
  fcs.tsne.cluster.figure.dir,
  fcs.umap.figure.dir,
  fcs.umap.cluster.figure.dir,
  trex.figure.dir,
  fcs.pca.figure.dir,
  fcs.mfi.stats.dir,
  fcs.cluster.stats.dir
)

table.dir <- fcs.cluster.table.dir

for ( the.dir in c( figure.dir, table.dir ) )
  if ( ! file.exists( the.dir ) )
    dir.create( the.dir, recursive = TRUE )


# function to set random seed depending on base number and string---------------

set.seed.here <- function( seed.base, seed.char )
{
  seed.add <- strtoi( substr( digest( seed.char, "xxhash32" ), 2, 8 ), 16 )
  seed.new <- seed.base + seed.add    
  set.seed( seed.new )
  invisible( seed.new )
}

# load data =========================

analysis.start.time <- Sys.time()

if (input.file.type == 2){
  file.extension <- "\\.fcs$"
}else {
  file.extension <- "\\.csv$"
}

flow.data.filename.all <- list.files( data.dir, file.extension )

flow.data.filename <- grep( paste0( fcs.condition, collapse = "|" ), 
                            flow.data.filename.all, value = TRUE )

sample.name.format <- paste0( "%s.%0", fcs.sample.number.width, "d" )
event.name.format <- paste0( "%s.%0", fcs.event.number.width, "d" )

flow.data.filename.sample <- rep( "", length( flow.data.filename ) )
names( flow.data.filename.sample ) <- flow.data.filename

sample.idx.next <- rep( 1, fcs.condition.n )
names( sample.idx.next ) <- fcs.condition

if (input.file.type==2){
  
  flow.data <- lapply( flow.data.filename, function( flow.data.fn ) {
    
    sample.flow.frame <- read.FCS( file.path( data.dir, flow.data.fn ), 
                                   transformation = NULL, truncate_max_range = FALSE )
    
    condition <- fcs.condition[ sapply( fcs.condition, grepl, flow.data.fn ) ]
    stopifnot( length( condition ) == 1 )
    
    sample.data <- exprs( sample.flow.frame )
    
    if ( ! all( fcs.channel %in% colnames( sample.data ) ) )
    {
      cat( sprintf( "File: %s\n", flow.data.fn ) )
      print( sort( fcs.channel[ 
        ! fcs.channel %in% colnames( sample.data ) ] ) )
      print( sort( colnames( sample.data ) ) )
      stop( "mismatch in names of fcs channels" )
    }
    
    sample.name <- sprintf( sample.name.format, condition, 
                            sample.idx.next[ condition ] )
    
    sample.data <- sample.data[ , fcs.channel, drop = FALSE ]
    
    event.n <- nrow( sample.data )
    if ( event.n > 0 ) {
      event.name <- sprintf( event.name.format, sample.name, 1 : event.n )
      rownames( sample.data ) <- event.name
    }
    
    flow.data.filename.sample[ flow.data.fn ] <<- sample.name
    sample.idx.next[ condition ] <<- sample.idx.next[ condition ] + 1
    
    sample.data
  } )
  
}else{
  flow.data <- lapply( flow.data.filename, function( flow.data.fn ) {
    
    sample.data <- as.matrix(fread( file.path( data.dir, flow.data.fn ), check.names = TRUE ))
    
    condition <- fcs.condition[ sapply( fcs.condition, grepl, flow.data.fn ) ]
    stopifnot( length( condition ) == 1 )
    
    if ( ! all( fcs.channel %in% colnames( sample.data ) ) )
    {
      cat( sprintf( "File: %s\n", flow.data.fn ) )
      print( sort( fcs.channel[ 
        ! fcs.channel %in% colnames( sample.data ) ] ) )
      print( sort( colnames( sample.data ) ) )
      cat("Channel mismatch error\n
        Please check that your files were all run with the same flow panel and try again.\n")
      Sys.sleep(message.delay.time*2)
      stop( "mismatch in names of channels" )
    }
    
    sample.name <- sprintf( sample.name.format, condition, 
                            sample.idx.next[ condition ] )
    
    sample.data <- sample.data[ , fcs.channel, drop = FALSE ]
    
    event.n <- nrow( sample.data )
    if ( event.n > 0 ) {
      event.name <- sprintf( event.name.format, sample.name, 1 : event.n )
      rownames( sample.data ) <- event.name
    }
    
    flow.data.filename.sample[ flow.data.fn ] <<- sample.name
    sample.idx.next[ condition ] <<- sample.idx.next[ condition ] + 1
    
    sample.data
  } )
}

flow.data <- do.call( rbind, flow.data )

if (input.file.type==2){
  colnames(flow.data) <- fcs.channel.label
}

# define samples---------------
flow.sample <- flow.data.filename.sample
names( flow.sample ) <- NULL

stopifnot( flow.sample == 
             unique( sub( "\\.[0-9]+$", "", rownames( flow.data ) ) ) )

flow.sample.n <- length ( flow.sample )

flow.sample.condition <- factor( sub( "\\.[0-9]+$", "", flow.sample ), 
                                 levels = fcs.condition )
names( flow.sample.condition ) <- flow.sample

# reorder samples to follow order of conditions
flow.sample <- flow.sample[ order( flow.sample.condition ) ]
flow.sample.condition <- flow.sample.condition[ flow.sample ]

flow.sample.label <- sapply( flow.sample, function( fs ) {
  sample.cond <- sub( "^(.*)\\.[0-9]+$", "\\1", fs )
  sample.num <- sub( "^.*\\.([0-9]+)$", "\\1", fs )
  sprintf( "%s-%s", fcs.condition.label[ sample.cond ], sample.num )
} )

flow.sample.filename <- sapply( flow.sample, function( fs  ) 
  names( which( flow.data.filename.sample == fs ) ) )

# define events
flow.event <- rownames( flow.data )
flow.event.n <- length( flow.event )

flow.event.sample <- factor( sub( "\\.[0-9]+$", "", flow.event ), 
                             levels = flow.sample )
names( flow.event.sample ) <- flow.event

flow.event.condition <- factor( sub( "\\.[0-9]+$", "", flow.event.sample ), 
                                levels = fcs.condition )
names( flow.event.condition ) <- flow.event

# reorder events to follow order of samples
flow.event.order <- order( flow.event.sample )

flow.data <- flow.data[ flow.event.order, ]
flow.event <- flow.event[ flow.event.order ]
flow.event.sample <- flow.event.sample[ flow.event.order ]
flow.event.condition <- flow.event.condition[ flow.event.order ]

flow.event.sample.n <- as.vector( table( flow.event.sample ) )
names( flow.event.sample.n ) <- flow.sample

flow.event.condition.n <- as.vector( table( flow.event.condition ) )
names( flow.event.condition.n ) <- fcs.condition.label

table( flow.sample.condition )
flow.event.condition.n
flow.event.sample.n

# define figure parameters for samples

flow.sample.color <- fcs.condition.color[ flow.sample.condition ]
names( flow.sample.color ) <- flow.sample

flow.sample.color.single <- unlist( lapply( fcs.condition, function( fc ) {
  cond.sample.n <- sum( flow.sample.condition == fc )
  rep( 
    fcs.color.pool, 
    ceiling( cond.sample.n / fcs.color.pool.n ) 
  )[ 1 : cond.sample.n ]
} ) )
names( flow.sample.color.single ) <- flow.sample

flow.sample.line.type <- fcs.condition.line.type[ flow.sample.condition ]
names( flow.sample.line.type ) <- flow.sample

flow.sample.line.type.single <- unlist( lapply( fcs.condition, function( fc ) {
  cond.sample.n <- sum( flow.sample.condition == fc )
  rep( 
    fcs.line.type.pool, 
    ceiling( cond.sample.n / fcs.line.type.pool.n ) 
  )[ 1 : cond.sample.n ]
} ) )
names( flow.sample.line.type.single ) <- flow.sample

flow.ce.diff.figure.dendrogram.weight.sample <- 
  fcs.ce.diff.figure.dendrogram.weight.condition[ flow.sample.condition ]
names( flow.ce.diff.figure.dendrogram.weight.sample ) <- flow.sample

# select data for dimensionality reduction---------------

set.seed.here( fcs.seed.base, "select data for dimensionality reduction" )

{
  if ( ! is.null( fcs.dmrd.data.sample.n ) )
  {
    if ( fcs.dmrd.data.sample.n < flow.event.n )
      dmrd.data.idx <- sort( sample( flow.event.n, fcs.dmrd.data.sample.n ) )
    else
      dmrd.data.idx <- 1 : flow.event.n
  }
  else if ( ! is.null( fcs.dmrd.data.sample.n.per.condition ) )
  {
    dmrd.data.idx <- unlist( sapply( fcs.condition, function( fc ) {
      fc.idx <- which( flow.event.condition == fc )
      if ( fcs.dmrd.data.sample.n.per.condition < length( fc.idx ) )
        sort( sample( fc.idx, fcs.dmrd.data.sample.n.per.condition ) )
      else
        fc.idx
    } ) )
    names( dmrd.data.idx ) <- NULL
  }
  else if ( ! is.null( fcs.dmrd.data.sample.n.per.sample ) )
  {
    dmrd.data.idx <- unlist( sapply( flow.sample, function( fs ) {
      fs.idx <- which( flow.event.sample == fs )
      if ( fcs.dmrd.data.sample.n.per.sample < length( fs.idx ) )
        sort( sample( fs.idx, fcs.dmrd.data.sample.n.per.sample ) )
      else
        fs.idx
    } ) )
    names( dmrd.data.idx ) <- NULL
  }
  else
    dmrd.data.idx <- 1 : flow.event.n
}

dmrd.data <- flow.data[ dmrd.data.idx, ]
dmrd.event.n <- dmrd.data.n <- nrow(dmrd.data)
dmrd.event.sample <- flow.event.sample[ dmrd.data.idx ]
dmrd.event.condition <- flow.event.condition[ dmrd.data.idx ]

# transform data if using fcs files----------
## TBD: include transforms for FACSDiscover, ZE5, Attune, Bigfoot
if ( input.file.type == 2 ){
  
  if( flowFrame@description$`$CYT` == "Aurora" ){
    width.basis <- -1000
    max.value <- 4194303
    log.decades <- 5.5
  } else if( flowFrame@description$`$CYT` == "ID7000" ){
    width.basis <- -500
    max.value <- 1000000
    log.decades <- 5
  } else {
    width.basis <- -100
    max.value <- 262144
    log.decades <- 4.5
  }
  
  extra.neg.decades <- 0
  
  biexp.transform <- flowjo_biexp(channelRange = 1250, maxValue = max.value, 
                                  pos = log.decades, neg = extra.neg.decades,
                                  widthBasis = width.basis)
  
  # transform data
  dmrd.data.untransformed <- dmrd.data
  dmrd.data <- apply( dmrd.data, 2, biexp.transform)
  dmrd.data <- apply( dmrd.data, 2, FUN = "-", 250 )
  dmrd.data.long <- pivot_longer(data.frame( dmrd.data ), cols = everything(), 
                                 names_to = "parameter", values_to = "value")
  rownames(dmrd.data) <- rownames(dmrd.data.untransformed)
  
  #plot histograms of channels
  transformation.plot <- ggplot( dmrd.data.long, 
                                 aes(x = value, y = after_stat(count) ))+
    geom_density(fill='black', alpha = 0.4) +
    theme_classic()+
    facet_wrap(~parameter, scales = "free")+
    coord_cartesian(xlim = c(-50,1000))+
    xlab("Channel")
  
  ggsave(
    file.path( fcs.density.figure.dir, 
               "biexponential_transform_fcs.jpg" ), 
    transformation.plot, 
    width = fcs.density.figure.width.base * ( fcs.channel.n + 1 )*1.2, 
    height = fcs.density.figure.height.base *100
  )
}

rm(flow.data)

## Clustering ===================

# calculate umap representation for data visualization---------------

{
  if ( fcs.use.cached.results && file.exists( fcs.umap.cache.file.path ) )
  {
    cat( "Using cached results for umap\n" )
    
    load( fcs.umap.cache.file.path )
  }
  else { 
    set.seed.here( fcs.seed.base, "calculate umap representation" )
    
    cat( "Calculating UMAP\n" )
    
    umap.result <- uwot::umap(dmrd.data, n_neighbors = fcs.tsne.perplexity,
                              n_epochs = fcs.umap.iter.n, 
                              n_threads = fcs.tsne.threads.n,
                              n_sgd_threads = fcs.tsne.threads.n, 
                              batch = TRUE, verbose = TRUE, ret_model = TRUE,
                              ret_nn = TRUE, ret_extra = "sigma")
    
    save(umap.result, file = fcs.umap.cache.file.path )
  }
}

umap.data <- umap.result$embedding
dimnames( umap.data ) <- NULL

dmrd.knn <- umap.result$nn$euclidean

umap.islands.plot <- ggplot(data.frame(umap.data), aes(x = X1, y = X2 )) + 
  geom_scattermore( alpha = fcs.dmrd.color.alpha, 
                    pointsize = fcs.umap.figure.point.size*2,
                    pixels = c(700,700))+
  theme_bw() + 
  theme( axis.title = element_blank(), 
         axis.text  = element_blank(), 
         axis.ticks = element_blank(), 
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         strip.background = element_blank(),
         plot.background = element_blank() )

ggsave( 
  file.path( fcs.umap.figure.dir, "umap_islands_plot.jpeg" ), 
  umap.islands.plot, 
  width = fcs.dmrd.figure.width, height = fcs.dmrd.figure.height 
)

# clustering with chosen method--------------

if (clustering.method == 2){
  # get FlowSOM clusters---------------
  # build som objects--embedsom method
  cat("Clustering data with FlowSOM via EmbedSOM\n")
  set.seed.here( fcs.seed.base, "get flowsom clusters" )
  
  flow.som <- EmbedSOM::SOM(dmrd.data, xdim = fcs.flow.som.dim, 
                            ydim = fcs.flow.som.dim, batch = TRUE,
                            parallel = TRUE, threads = fcs.tsne.threads.n )
  
  # get clusters
  flow.som.mapping <- flow.som$mapping[ , 1 ]
  flow.som.codes <- flow.som$codes
  
  # get clusters from som mapping
  consensus.cluster <- ConsensusClusterPlus( t( flow.som.codes ),
                                             maxK = fcs.cluster.n, reps = 100, pItem = 0.9, pFeature = 1,
                                             clusterAlg = "hc", innerLinkage = "average", finalLinkage = "average",
                                             distance = "euclidean", 
                                             seed = set.seed.here( fcs.seed.base, "get clusters from som mapping" ) )
  
  flow.som.event.cluster <- consensus.cluster[[ fcs.cluster.n ]]$
    consensusClass[ flow.som.mapping ]
  
  # reorder clusters from bigger to smaller
  flow.som.cluster.rank <- 1 + fcs.cluster.n - 
    rank( table( flow.som.event.cluster ), ties.method = "last" )
  flow.som.event.cluster <- flow.som.cluster.rank[ flow.som.event.cluster ]
  names( flow.som.event.cluster ) <- NULL
  
  # set clusters as a factor
  flow.som.event.cluster <- factor( flow.som.event.cluster, 
                                    levels = 1 : fcs.cluster.n )
  levels( flow.som.event.cluster ) <- fcs.cluster
  
  dmrd.event.cluster <- flow.som.event.cluster
  dmrd.event.cluster.n <- as.vector( table( dmrd.event.cluster ) )
  names( dmrd.event.cluster.n ) <- fcs.cluster
  
  length( dmrd.event.cluster )
  table( dmrd.event.cluster )
  
} else {
  # get phenograph clusters---------------
  cat("Clustering data with Phenograph\n")
  
  set.seed.here( fcs.seed.base, "get clusters with Phenograph" )
  links <- FastPG::rcpp_parallel_jce( dmrd.knn$idx )
  links <- FastPG::dedup_links( links )
  clusters <- FastPG::parallel_louvain( links )
  
  phenograph.event.cluster <- factor(clusters$communities)
  
  fcs.cluster.n <- length(unique(phenograph.event.cluster))
  phenograph.cluster.rank <- 1 + fcs.cluster.n - 
    rank( table( phenograph.event.cluster ), ties.method = "last" )
  phenograph.event.cluster <- phenograph.cluster.rank[ phenograph.event.cluster ]
  names( phenograph.event.cluster ) <- NULL
  
  # set clusters as a factor
  phenograph.event.cluster <- factor( phenograph.event.cluster, 
                                      levels = 1 : fcs.cluster.n )
  levels( phenograph.event.cluster ) <- fcs.cluster
  
  # reorder events
  dmrd.event.cluster <- phenograph.event.cluster
  
  dmrd.event.cluster.n <- as.vector( table( dmrd.event.cluster ) )
  names( dmrd.event.cluster.n ) <- fcs.cluster
  
  length( dmrd.event.cluster )
  table( dmrd.event.cluster )
}

## Automated cluster naming ===============

# select correct databases based on species
cell.marker.filename <- paste0(species.used, "_marker_names.xlsx")
cell.type.filename <- paste0(species.used, "_celltype_database.xlsx")

cell.database <- read_xlsx( file.path(src.dir, cell.type.filename) )

# restrict cell types based on tissues used
tissue.type <- unique(cell.database$Tissue.restricted)[tissue.type]

if ( length(tissue.type) > 1 ){
  cell.database <- dplyr::filter(cell.database, Tissue.restricted %in% tissue.type)
} else if ( tissue.type == 0 ){
  cell.database <- cell.database
} else {
  cell.database <- dplyr::filter(cell.database, Tissue.restricted %in% tissue.type)
}

# restrict cell ID based on known input
selected.cell.type <- unique(cell.database$Parent.cell.type)[selected.cell.type]


# match markers to database
marker.synonyms <- read_xlsx( file.path(src.dir, cell.marker.filename) )

marker.alternates <- lapply(1:nrow(marker.synonyms), 
                            function(x) gsub(" ","",unlist(strsplit(toString(marker.synonyms$Alternates[x]),","))))
marker.alternates <- lapply(1:nrow(marker.synonyms), function(x)
  gsub("\\*"," ", unlist(marker.alternates[x])))
names(marker.alternates) <- marker.synonyms$Marker.Name

new.marker.names <- lapply( 1:length(fcs.channel.label), function(x){
  names(marker.alternates)[grep( paste0( fcs.channel.label[x], "\\b"), marker.alternates,
                                 ignore.case = TRUE)]
} )

# scale data for cluster matching
flow.data.cluster.median <- apply( dmrd.data, 2, tapply, 
                                   dmrd.event.cluster, median )

scaled.fcs.cluster.median <- scale(flow.data.cluster.median, scale = FALSE)

colnames(scaled.fcs.cluster.median) <- new.marker.names

# prepare marker lists and match clusters
source( file.path( src.dir, "prepare_marker_lists.r" ) )
source( file.path( src.dir, "flow_cluster_id_score.r" ) )

marker.list <- prepare_marker_lists( file.path( src.dir, cell.type.filename ), 
                                     tissue.type, selected.cell.type )

auto.fcs.cluster.label <- fcs.cluster.label

# match clusters to cell type database
scaled.id.score <- flow_cluster_id_score(t(scaled.fcs.cluster.median),
                                         marker_pos = marker.list$markers_positive, 
                                         marker_neg = marker.list$markers_negative )

for (cluster in 1:ncol(scaled.id.score)){
  auto.fcs.cluster.label[cluster] <- names( which.max(scaled.id.score[,cluster]) )
}

# plot heatmaps of cluster ID scores--------

jpeg( file.path(fcs.density.figure.dir, "cluster_id_heatmap.jpg"), 
      width = 1000, height = 1000 )
heatmap(scaled.id.score, Rowv = NA, Colv = NA, scale = "none",
        margins = c(5,10),
        xlab = "Clusters", ylab = "matching cell types")
dev.off()

jpeg( file.path(fcs.density.figure.dir, "cluster_id_heatmap_dendro.jpg"), 
      width = 1000, height = 1000 )
heatmap(scaled.id.score, scale = "none",
        margins = c(5,10),
        xlab = "Clusters", ylab = "matching cell types")
dev.off()



# differentiate similar clusters by variable markers

rownames(flow.data.cluster.median) <- auto.fcs.cluster.label
colnames(flow.data.cluster.median) <- fcs.channel.label

# find (positions of) clusters with non-unique names
non.unique.cluster.names <- auto.fcs.cluster.label[duplicated(auto.fcs.cluster.label)]
# sort into groups
non.unique.groups <- unique(non.unique.cluster.names)

# append high variance marker labels to all clusters until unique

for (cluster.group in non.unique.groups) {
  
  # collect clusters needing renaming
  data.temp <- flow.data.cluster.median[ grepl( cluster.group, rownames(flow.data.cluster.median)),]
  
  # find most variable channels
  variances <- apply(data.temp, 2, var)
  sorted.variances <- sort(variances, decreasing = TRUE)
  
  # find min and max (neg and pos) for each channel
  max.expression <- apply(data.temp, 2, max)
  min.expression <- apply(data.temp, 2, min)
  
  for( cluster in 1:nrow(data.temp)){
    # get position of cluster in full list
    position.in.cluster.list <- rownames(flow.data.cluster.median)[
      grepl(cluster.group, rownames(flow.data.cluster.median))][cluster]
    
    # append markers to name for each positive variable up to the number of clusters in the group or max 4
    
    markers.to.append.n <- ifelse( nrow(data.temp) < 5, nrow(data.temp), 4 )
    
    for( marker in names(sorted.variances)[1:markers.to.append.n]){
      rownames(flow.data.cluster.median)[as.numeric(names(position.in.cluster.list))] <- 
        if ( abs(flow.data.cluster.median[as.numeric(names(position.in.cluster.list)),marker] - max.expression[marker] ) < abs(flow.data.cluster.median[as.numeric(names(position.in.cluster.list)),marker] - min.expression[marker] )){
          paste( rownames(flow.data.cluster.median)[as.numeric(names(position.in.cluster.list))], marker )
        } else{
          rownames(flow.data.cluster.median)[as.numeric(names(position.in.cluster.list))] <- rownames(flow.data.cluster.median)[as.numeric(names(position.in.cluster.list))]
        }
    }
    
  }
  
}

fcs.cluster.label <- rownames(flow.data.cluster.median)

# plot histograms of marker expression

source( file.path( src.dir, "plot_data_histograms.r") )

source( file.path( src.dir, "plot_cluster_histograms.r") )


## Optional: rename clusters in this section============
# run command below and paste output into flowcytoscript_parameter.r on line 141
# where it also says fcs.cluster.label
# change the names to suit you in the parameter file, not here

cat('c("', paste0(fcs.cluster.label, collapse = '","'), '")', sep = "")

param.filename <- "./flowcytoscript_parameter.r"
source( param.filename )

# run command below and paste output into flowcytoscript_parameter.r on line 141
# where it also says fcs.cluster.color
# change the colors to suit you in the parameter file
# check the R Color Guide for color names
cat('c("', paste0(fcs.cluster.color, collapse = '","'), '")', sep = "")

param.filename <- "./flowcytoscript_parameter.r"
source( param.filename )

# plot cluster histograms with new colors and names
source( file.path( src.dir, "plot_cluster_histograms.r") )


## end manual cluster renaming===============




# save cluster counts

flow.cluster.count <- sapply( fcs.cluster, function( fc )
  table( dmrd.event.sample[ dmrd.event.cluster == fc ] ) )

stopifnot( rownames( flow.cluster.count ) == flow.sample )
stopifnot( colnames( flow.cluster.count ) == fcs.cluster )

rownames( flow.cluster.count ) <- flow.sample.label
colnames( flow.cluster.count ) <- fcs.cluster.label

write.csv( flow.cluster.count, file = file.path( fcs.cluster.table.dir, 
                                                 sprintf( "%s.csv", fcs.cluster.table.counts ) ) )

flow.cluster.percent <- flow.cluster.count/rowSums(flow.cluster.count)*100

write.csv( flow.cluster.percent, file = file.path( fcs.cluster.table.dir, 
                                                   sprintf( "%s.csv", "cluster_proportions" ) ) )


# plot heatmaps---------------

heatmap.type <- c( "by_condition", "by_sample", "by_cluster" )

for ( ht in heatmap.type )
{
  if ( ht == "by_condition" ) {
    flow.data.group.median <- apply( dmrd.data, 2, tapply, 
                                     dmrd.event.condition, median )
    margin.col <- fcs.heatmap.label.factor.col * 
      max( nchar( fcs.condition.label ) )
    group.label <- fcs.condition.label
    group.color <- fcs.condition.color
  }
  else if ( ht == "by_sample" ) {
    flow.data.group.median <- apply( dmrd.data, 2, tapply, 
                                     dmrd.event.sample, median )
    margin.col <- fcs.heatmap.label.factor.col * 
      max( nchar( flow.sample.label ) )
    group.label <- flow.sample.label
    group.color <- flow.sample.color
  }
  else if ( ht == "by_cluster" ) {
    flow.data.group.median <- apply( dmrd.data, 2, tapply, 
                                     dmrd.event.cluster, median )
    margin.col <- fcs.heatmap.label.factor.col * 
      max( nchar( fcs.cluster.label ) )
    group.label <- fcs.cluster.label
    group.color <- fcs.cluster.color
  }
  else
    stop( "wrong heatmap type" )
  
  margin.row <- 
    fcs.heatmap.label.factor.row * max( nchar( fcs.channel.label ) )
  
  if ( margin.row < 5 )
    margin.row <- 5
  if ( margin.col < 5 )
    margin.col <- 5
  
  png( filename = file.path( fcs.heatmap.figure.dir, 
                             sprintf( "%s_%s.png", fcs.heatmap.figure, ht ) ), 
       width = fcs.heatmap.width, height = fcs.heatmap.height )
  heatmap( t( flow.data.group.median ), scale = "row", 
           labRow = fcs.channel.label,  labCol = group.label, 
           col = fcs.heatmap.palette, ColSideColors = group.color, 
           cexRow = fcs.heatmap.font.size, cexCol = fcs.heatmap.font.size, 
           margins = c( margin.col, margin.row ) )
  dev.off()
}


# plot histograms---------------

histogram.type <- c( "by_cluster", "by_sample", "by_assay" )

condition.sample.n <- as.vector( table( flow.sample.condition ) )

dmrd.event.sample.n <- as.numeric( table( dmrd.event.sample ) )

for ( ht in histogram.type )
{
  flow.data.cluster.fraction <- lapply( fcs.cluster, function( fc ) {
    sample.event.n <- as.vector( table( 
      dmrd.event.sample[ dmrd.event.cluster == fc ] ) )
    
    if ( ht == "by_cluster" )
    {
      sample.fraction <- log2( 
        ( sample.event.n / dmrd.event.cluster.n[ fc ] ) / 
          ( dmrd.event.sample.n / dmrd.data.n ) 
      )
      sample.fraction[ is.infinite( sample.fraction ) ] <- NA
      fraction <- tapply( sample.fraction, flow.sample.condition, mean, 
                          na.rm = TRUE )
      std.err <- tapply( sample.fraction, flow.sample.condition, sd, 
                         na.rm = TRUE ) / sqrt( condition.sample.n )
    }
    else if ( ht == "by_sample" )
    {
      sample.fraction <- 100 * sample.event.n / dmrd.event.sample.n
      fraction <- tapply( sample.fraction, flow.sample.condition, mean, 
                          na.rm = TRUE )
      std.err <- tapply( sample.fraction, flow.sample.condition, sd, 
                         na.rm = TRUE ) / sqrt( condition.sample.n )
    }
    else if ( ht == "by_assay" )
    {
      sample.fraction <- 100 * sample.event.n / dmrd.data.n
      fraction <- tapply( sample.fraction, flow.sample.condition, sum, 
                          na.rm = TRUE )
      std.err <- NA
    }
    else
      stop( "wrong histogram type" )
    
    data.frame( cluster = fc, condition = fcs.condition, fraction, std.err )
  } )
  
  flow.data.cluster.fraction <- do.call( rbind, flow.data.cluster.fraction )
  flow.data.cluster.fraction$condition <- factor( 
    flow.data.cluster.fraction$condition, levels = fcs.condition )
  
  histogram.plot <- ggplot( flow.data.cluster.fraction, 
                            aes( x = cluster, y = fraction, fill = condition ) ) + 
    geom_bar( stat = "identity", position = position_dodge2() ) + 
    geom_errorbar( aes( ymin = fraction - std.err, 
                        ymax = fraction + std.err ), 
                   linewidth = fcs.histogram.error.bar.size,
                   position = position_dodge2() ) + 
    scale_x_discrete( labels = fcs.cluster.label, name = "" ) + 
    scale_fill_manual( values = fcs.condition.color, 
                       limits = fcs.condition, breaks = fcs.condition, 
                       labels = fcs.condition.label ) + 
    theme_bw() + 
    theme( 
      axis.text = element_text( size = fcs.histogram.font.size, 
                                angle = 90 ), 
      axis.title = element_text( size = fcs.histogram.font.size + 1 ), 
      legend.title = element_blank(), 
      legend.text = element_text( size = fcs.histogram.font.size + 1 ), 
      legend.key.size = unit( fcs.histogram.legend.key.size, "lines" ), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank() 
    )
  
  if ( ht == "by_cluster" )
    histogram.plot <- histogram.plot + 
    ylab( "Average log2( frequency ratio )" )
  else if ( ht == "by_sample" )
    histogram.plot <- histogram.plot + ylab( "Average frequency" )
  else if ( ht == "by_assay" )
    histogram.plot <- histogram.plot + ylab( "Frequency" )
  
  ggsave( 
    file.path( fcs.histogram.figure.dir, 
               sprintf( "%s_%s.png", fcs.histogram.figure, ht ) ), 
    histogram.plot, 
    width = fcs.histogram.width, 
    height = fcs.histogram.height + 
      fcs.histogram.label.factor.height * 
      max( nchar( fcs.cluster.label ) ) 
  )
}


# run pca on median values for each marker for each sample------------

flow.data.group.median <- data.frame(apply( dmrd.data, 2, tapply, 
                                            dmrd.event.sample, median ))

pca.mfi.result <- prcomp(flow.data.group.median)

# add labels and extract coordinates
flow.data.group.median$Sample <- rownames(flow.data.group.median)
flow.data.group.median$sample.color <- flow.sample.color[flow.data.group.median$Sample]
flow.data.group.median$Group <- flow.sample.condition[flow.data.group.median$Sample]
flow.data.group.median[ , c("PC1", "PC2")] <- pca.mfi.result$x[, 1:2]
pca.mfi.loadings <- data.frame(pca.mfi.result$rotation)
pca.mfi.loadings$var <- rownames(pca.mfi.result$rotation)
mfi.scaling.loadings <- max(abs( flow.data.group.median[,c("PC1", "PC2")] )) / max( abs(pca.mfi.loadings[, 1:2])) / 2
pca.mfi.loadings[, 1:2] <- pca.mfi.loadings[, 1:2] * mfi.scaling.loadings
pc.mfi.variance <- summary(pca.mfi.result)
pc.mfi.variance <- data.frame(pc.mfi.variance$importance)
pc1.mfi.proportion <- pc.mfi.variance[2,1]*100
pc2.mfi.proportion <- pc.mfi.variance[2,2]*100

# plot pca on markers with loadings
pca.marker.var.plot <- ggplot(data = flow.data.group.median, 
                              aes(PC1, PC2, col = Group )) +
  geom_point( size = pca.point.size ) +
  scale_color_manual( values = fcs.condition.color,
                      breaks = fcs.condition,
                      labels = fcs.condition.label)+
  geom_segment( inherit.aes = FALSE, data = pca.mfi.loadings, 
                aes(x = 0, y = 0, xend = PC1, yend = PC2 ),
                color = "black", linewidth = pca.loading.arrow.size, 
                arrow = arrow(length = unit(0.03, "npc")) ) +
  geom_label( inherit.aes = FALSE, data = pca.mfi.loadings, 
              aes(PC1 * 1.2, PC2 * 1.2, label = pca.mfi.loadings$var ), 
              size = pca.loading.text.size,
              label.size = pca.loading.label.size ) +
  xlab(paste( "PC1 (", pc1.mfi.proportion, "%)", sep = "" )) +
  ylab(paste( "PC2 (", pc2.mfi.proportion, "%)", sep = "" )) +
  theme_bw() +
  theme( panel.grid = element_blank(), legend.title = element_blank() )

ggsave(filename = file.path( fcs.pca.figure.dir, 
                             paste0( pca.figure.mfi, "_loadings.png" )),
       pca.marker.var.plot, width = fcs.pca.figure.width, height = fcs.pca.figure.height )

# plot pca on markers without loadings
pca.marker.sample.plot <- ggplot(data = flow.data.group.median, 
                                 aes(PC1, PC2, col = Group)) +
  geom_point(size = pca.point.size ) +
  scale_color_manual( values = fcs.condition.color,
                      breaks = fcs.condition,
                      labels = fcs.condition.label)+
  xlab(paste( "PC1 (", pc1.mfi.proportion, "%)", sep = "" )) +
  ylab(paste( "PC2 (", pc2.mfi.proportion, "%)", sep = "" ))+
  theme_bw() +
  theme(panel.grid = element_blank(), legend.title = element_blank())

ggsave(filename = file.path( fcs.pca.figure.dir, 
                             paste( pca.figure.mfi, "_samples.png" )),
       pca.marker.sample.plot, width = fcs.pca.figure.width, height = fcs.pca.figure.height )

# calculate inter-group distances for marker PCA-------------

group.mfi.pca.data <- data.frame(pca.mfi.result$x)
group.mfi.pca.data$Group <- flow.sample.condition[flow.data.group.median$Sample]
group.mfi.pca.data$Group <- fcs.condition.label[group.mfi.pca.data$Group]

group.mfi.summary <- group.mfi.pca.data %>% 
  dplyr::group_by(Group) %>%
  dplyr::summarise_all( median )

group.mfi.distance <- dist(group.mfi.summary[,-1])
group.mfi.distance <- as.matrix(group.mfi.distance, labels=TRUE)
group.mfi.distance <- as.data.frame(group.mfi.distance)
colnames(group.mfi.distance) <- rownames(group.mfi.distance) <- group.mfi.summary[['Group']]

write.csv(group.mfi.distance, file = file.path(fcs.pca.figure.dir, "intergroup_mfi_distance.csv"))

png( filename = file.path( fcs.pca.figure.dir, 
                           "pca_marker_distance_dendrogram.png" ), 
     width = fcs.heatmap.width/2, height = fcs.heatmap.height/2 )
heatmap(as.matrix(group.mfi.distance),
        col = fcs.heatmap.palette,
        margins = c( margin.col/2, margin.row ) )
legend(x = "topright", legend = c("close", "far"),
       fill = fcs.heatmap.palette[c(1,100)], cex = 2 )
dev.off()

# run pca on cluster frequencies------------
cat("Running PCA on cluster distributions\n")

pca.cluster.result <- prcomp(flow.cluster.percent)

# add labels and extract coordinates

flow.cluster.frame <- data.frame(flow.cluster.percent)

flow.cluster.frame$Sample <- names(flow.sample.color)

flow.cluster.frame$sample.color <- flow.sample.color[flow.cluster.frame$Sample]

flow.cluster.frame$Group <- flow.data.group.median$Group

flow.cluster.frame[ , c("PC1", "PC2")] <- pca.cluster.result$x[, 1:2]
pca.cluster.loadings <- data.frame(pca.cluster.result$rotation)
pca.cluster.loadings$var <- rownames(pca.cluster.result$rotation)
cluster.scaling.loadings <- max(abs( flow.cluster.frame[,c("PC1", "PC2")] )) / max( abs(pca.cluster.loadings[, 1:2])) / 2
pca.cluster.loadings[, 1:2] <- pca.cluster.loadings[, 1:2] * cluster.scaling.loadings
pc.cluster.variance <- summary(pca.cluster.result)
pc.cluster.variance <- data.frame(pc.cluster.variance$importance)
pc1.cluster.proportion <- pc.cluster.variance[2,1]*100
pc2.cluster.proportion <- pc.cluster.variance[2,2]*100

# plot pca on clusters with loadings
pca.cluster.var.plot <- ggplot(data = flow.cluster.frame, 
                               aes(PC1, PC2, col = Group )) +
  geom_point( size = pca.point.size ) +
  scale_color_manual( values = fcs.condition.color,
                      breaks = fcs.condition,
                      labels = fcs.condition.label)+
  geom_segment( inherit.aes = FALSE, data = pca.cluster.loadings, 
                aes(x = 0, y = 0, xend = PC1, yend = PC2 ),
                color = "black", linewidth = pca.loading.arrow.size, 
                arrow = arrow(length = unit(0.03, "npc")) ) +
  geom_label( inherit.aes = FALSE, data = pca.cluster.loadings, 
              aes(PC1 * 1.2, PC2 * 1.2, label = pca.cluster.loadings$var ), 
              size = pca.loading.text.size,
              label.size = pca.loading.label.size ) +
  xlab(paste( "PC1 (", pc1.cluster.proportion, "%)", sep = "" )) +
  ylab(paste( "PC2 (", pc2.cluster.proportion, "%)", sep = "" )) +
  theme_bw() +
  theme( panel.grid = element_blank(), legend.title = element_blank() )

ggsave(filename = file.path( fcs.pca.figure.dir, 
                             paste0( pca.figure.cluster, "_loadings.png" )),
       pca.cluster.var.plot, width = fcs.pca.figure.width, height = fcs.pca.figure.height )


# plot pca on clusters without loadings
pca.cluster.sample.plot <- ggplot(data = flow.cluster.frame, 
                                  aes(PC1, PC2, col = Group )) +
  geom_point(size = pca.point.size ) +
  scale_color_manual( values = fcs.condition.color,
                      breaks = fcs.condition,
                      labels = fcs.condition.label)+
  xlab(paste( "PC1 (", pc1.cluster.proportion, "%)", sep = "" )) +
  ylab(paste( "PC2 (", pc2.cluster.proportion, "%)", sep = "" ))+
  theme_bw() +
  theme(panel.grid = element_blank(), legend.title = element_blank())

ggsave(filename = file.path( fcs.pca.figure.dir, 
                             paste( pca.figure.cluster, "_samples.png" )),
       pca.cluster.sample.plot, width = fcs.pca.figure.width, height = fcs.pca.figure.height )

# calculate inter-group distances for cluster PCA-------
group.cluster.pca.data <- data.frame(pca.cluster.result$x)
group.cluster.pca.data$Group <- flow.sample.condition[flow.data.group.median$Sample]
group.cluster.pca.data$Group <- fcs.condition.label[group.cluster.pca.data$Group]

group.cluster.summary <- group.cluster.pca.data %>% 
  dplyr::group_by(Group) %>%
  dplyr::summarise_all( median )

group.cluster.distance <- dist(group.cluster.summary[,-1])
group.cluster.distance <- as.matrix(group.cluster.distance, labels=TRUE)
group.cluster.distance <- as.data.frame(group.cluster.distance)
colnames(group.cluster.distance) <- rownames(group.cluster.distance) <- group.cluster.summary[['Group']]

write.csv(group.cluster.distance, file = file.path(fcs.pca.figure.dir, "intergroup_cluster_distance.csv"))

png( filename = file.path( fcs.pca.figure.dir, 
                           "pca_cluster_distance_dendrogram.png" ), 
     width = fcs.heatmap.width/2, height = fcs.heatmap.height/2 )
heatmap(as.matrix(group.cluster.distance),
        col = fcs.heatmap.palette,
        margins = c( margin.col/2, margin.row ) )
legend(x = "topright", legend = c("close", "far"),
       fill = fcs.heatmap.palette[c(1,100)], cex = 2 )
dev.off()


# run ANOVA with Tukey's post-test on marker MFIs and cluster frequencies-------

source( file.path( src.dir, "flowcytoscript_anova.r") )


# calculate tsne representation---------------

{
  if ( fcs.use.cached.results && file.exists( fcs.tsne.cache.file.path ) )
  {
    cat( "Using cached results for tsne\n" )
    
    load( fcs.tsne.cache.file.path )
  }
  else
  {
    fcs.tsne.learning.rate <- ifelse( nrow(dmrd.data)/fcs.tsne.exaggeration.factor > 2000, 
                                      nrow(dmrd.data)/fcs.tsne.exaggeration.factor, 2000 )
    
    if ( is.null(dmrd.knn) ){
      set.seed.here( fcs.seed.base, "calculate tsne representation" )
      
      cat( "Calculating tSNE\n" )
      
      tsne.result <- Rtsne( dmrd.data, perplexity = fcs.tsne.perplexity, 
                            exaggeration_factor = fcs.tsne.exaggeration.factor,
                            eta = fcs.tsne.learning.rate,
                            max_iter = fcs.tsne.iter.n, stop_lying_iter = fcs.tsne.early.iter.n,
                            check_duplicates = FALSE, pca = FALSE, 
                            num_threads = fcs.tsne.threads.n )
      
    } else {
      set.seed.here( fcs.seed.base, "calculate tsne representation" )
      
      cat( "Calculating tSNE\n" )
      
      tsne.result <- Rtsne_neighbors( dmrd.knn$idx, dmrd.knn$dist, 
                                      perplexity = fcs.tsne.perplexity, 
                                      exaggeration_factor = fcs.tsne.exaggeration.factor,
                                      eta = fcs.tsne.learning.rate,
                                      max_iter = fcs.tsne.iter.n, stop_lying_iter = fcs.tsne.early.iter.n,
                                      check_duplicates = FALSE, pca = FALSE, 
                                      num_threads = fcs.tsne.threads.n )
    }
    
    
    save( tsne.result, file = fcs.tsne.cache.file.path )
  }
}

tsne.data <- tsne.result$Y

# plot data =================================
# plot tsne convergence---------------

tsne.iter <- 1 + 50 * ( 1 : length( tsne.result$itercosts ) - 1 )
tsne.cost <- tsne.result$itercosts

png( filename = file.path( fcs.tsne.figure.dir, 
                           sprintf( "%s.png", fcs.tsne.figure.convergence ) ), 
     width = fcs.tsne.figure.convergence.width, 
     height = fcs.tsne.figure.convergence.height )
par( mar = c( 5, 5.8, 2, 1.4 ) )
plot( tsne.iter, tsne.cost, log = "y", xlab = "Iteration", 
      ylab = "Total cost", xlim = c( 0, fcs.tsne.iter.n ), 
      ylim = c( 1, max( tsne.cost ) ), col = "blue3", 
      pch = 20, cex = 1, cex.lab = 2.5, cex.axis = 2 )
abline( h = tsne.cost[ length( tsne.cost ) ], lty = 2, lwd = 1.5, 
        col = "blue3" )
dev.off()


# plot tsne representations---------------

tsne.data.max <- max( abs( tsne.data ) )

# plot tSNE just with clusters colored
if ( use.scattermore == TRUE ){
  plot.cluster.drmd.figures(
    tsne.data, tsne.data.max, 
    fcs.tsne.figure.lims.factor, fcs.tsne.figure.point.size, 
    fcs.tsne.cluster.figure.dir, fcs.tsne.figure.plot, 
    dmrd.data, dmrd.event.cluster, dmrd.event.condition )
} else {
  ggplot.cluster.drmd.figures(
    tsne.data, tsne.data.max, 
    fcs.tsne.figure.lims.factor, fcs.tsne.figure.point.size, 
    fcs.tsne.cluster.figure.dir, fcs.tsne.figure.plot, 
    dmrd.data, dmrd.event.cluster, dmrd.event.condition )
}


# plot tSNE with everything
if ( use.scattermore == TRUE ){
  plot.all.dmrd.figures( 
    tsne.data, tsne.data.max, 
    fcs.tsne.figure.lims.factor, fcs.tsne.figure.point.size, 
    fcs.tsne.figure.dir, fcs.tsne.figure.plot, 
    dmrd.data, dmrd.event.cluster, dmrd.event.condition )
} else {
  ggplot.all.dmrd.figures( 
    tsne.data, tsne.data.max, 
    fcs.tsne.figure.lims.factor, fcs.tsne.figure.point.size, 
    fcs.tsne.figure.dir, fcs.tsne.figure.plot, 
    dmrd.data, dmrd.event.cluster, dmrd.event.condition )
}


# plot umap representations---------------

umap.data.max <- max( abs(umap.data) )

# plot umap just with clusters colored
if ( use.scattermore == TRUE ){
  plot.cluster.drmd.figures(
    umap.data, umap.data.max, 
    fcs.umap.figure.lims.factor, fcs.umap.figure.point.size, 
    fcs.umap.cluster.figure.dir, fcs.umap.figure.plot, 
    dmrd.data, dmrd.event.cluster, dmrd.event.condition )
} else {
  ggplot.cluster.drmd.figures(
    umap.data, umap.data.max, 
    fcs.umap.figure.lims.factor, fcs.umap.figure.point.size, 
    fcs.umap.cluster.figure.dir, fcs.umap.figure.plot, 
    dmrd.data, dmrd.event.cluster, dmrd.event.condition )
}


# plot umap with everything
if ( use.scattermore == TRUE ){
  plot.all.dmrd.figures( 
    umap.data, umap.data.max, 
    fcs.umap.figure.lims.factor, fcs.umap.figure.point.size, 
    fcs.umap.figure.dir, fcs.umap.figure.plot, 
    dmrd.data, dmrd.event.cluster, dmrd.event.condition )
} else {
  ggplot.all.dmrd.figures( 
    umap.data, umap.data.max, 
    fcs.umap.figure.lims.factor, fcs.umap.figure.point.size, 
    fcs.umap.figure.dir, fcs.umap.figure.plot, 
    dmrd.data, dmrd.event.cluster, dmrd.event.condition )
  
}


# trex: highlight changed regions ================

# subset data based on trex groups
trex.selection <- match(trex.condition, fcs.condition)
trex.dmrd.data <- subset(dmrd.data, dmrd.event.condition %in% trex.condition)
trex.umap.data <- umap.data
rownames(trex.umap.data) <- rownames(dmrd.data)
trex.umap.data <- subset(trex.umap.data, dmrd.event.condition %in% trex.condition)
trex.tsne.data <- tsne.data
rownames(trex.tsne.data) <- rownames(dmrd.data)
trex.tsne.data <- subset(trex.tsne.data, dmrd.event.condition %in% trex.condition)

# calculate knn for the trex groups

if ( length(fcs.condition) != length(trex.condition) ) {
  set.seed.here( fcs.seed.base, "get knn with Rcpp" )
  trex.knn <- RcppHNSW::hnsw_knn( trex.dmrd.data, k= fcs.tsne.perplexity,
                                  distance= 'l2', n_threads= fcs.tsne.threads.n )
} else if ( !is.null(dmrd.knn) ) {
  trex.knn <- dmrd.knn
} else {
  set.seed.here( fcs.seed.base, "get knn with Rcpp" )
  trex.knn <- RcppHNSW::hnsw_knn( trex.dmrd.data, k= fcs.tsne.perplexity,
                                  distance= 'l2', n_threads= fcs.tsne.threads.n )
}

trex.knn.index <- trex.knn$idx
first.condition.length <- nrow(subset(dmrd.data, dmrd.event.condition == trex.condition[1]))
trex.knn.index[trex.knn.index <= first.condition.length] <- 0
trex.knn.index[trex.knn.index > first.condition.length] <- 1

# calculate percent change in each KNN region
percent.change.knn <- (rowSums(trex.knn.index) / fcs.tsne.perplexity * 100)


# create trex plots

plot.changed.regions(
  trex.tsne.data, percent.change.knn,
  trex.figure.point.size, trex.condition, 
  trex.figure.dir, fcs.tsne.figure.plot
)

plot.changed.regions(
  trex.umap.data, percent.change.knn,
  trex.figure.point.size, trex.condition,
  trex.figure.dir, fcs.umap.figure.plot
)

# calculate stats for T-REX groups

trex.group.cluster.long <- flow.cluster.frame.long %>%
  dplyr::filter( Group %in% trex.condition)

p.values.all <- NULL
p.values.all <- list()

for (cluster in unique(trex.group.cluster.long$Cluster)) {
  temp <- trex.group.cluster.long %>% 
    dplyr::filter(Cluster == cluster) %>%
    group_by(Group) %>%
    suppressMessages(summarise(Percent = Percent))
  
  res_aov <- aov( Percent ~ Group, data = temp )
  
  fitted.em <- emmeans(res_aov, "Group", data = temp )
  
  p.values <- data.frame( pairs(fitted.em, adjust = "tukey" ) )
  p.values <- p.values[,c(1,6)]
  p.values$significant <- ifelse( p.values$p.value < 0.05, "Yes", "No" )
  p.values.all[cluster] <- p.values$p.value
  
  write.csv(p.values, file.path( trex.figure.dir, paste0( cluster, "_frequencies_anova.csv")))
  
}


# run crossentropy if desired-----------

if(run.crossentropy == TRUE){
  # calculate cross-entropy test for tsne by condition---------------
  
  set.seed.here( fcs.seed.base, "calculate cross-entropy test for tsne" )
  
  ce.diff.test.tsne.res <- ce.diff.test.tsne( 
    dmrd.data, tsne.data, 
    dmrd.event.condition, 
    partition.label = fcs.condition.label, 
    partition.color = fcs.condition.color, 
    partition.line.type = fcs.condition.line.type, 
    base.test = fcs.ce.diff.base.test, 
    base.dist = fcs.ce.diff.base.dist, 
    prob.sample.n = fcs.ce.diff.prob.sample.n, 
    dendrogram.order.weight = fcs.ce.diff.figure.dendrogram.weight.condition, 
    result = file.path( fcs.ce.diff.tsne.figure.dir, 
                        sprintf( "%s_condition.txt", fcs.ce.diff.tsne.result ) ), 
    cdf.figure = file.path( fcs.ce.diff.tsne.figure.dir, 
                            sprintf( "%s_condition.png", fcs.ce.diff.tsne.figure.cdf ) ), 
    dendrogram.figure = file.path( fcs.ce.diff.tsne.figure.dir, 
                                   sprintf( "%s_condition.png", fcs.ce.diff.tsne.figure.dendrogram ) ) )
  
  
  # calculate cross-entropy test for tsne by sample---------------
  
  set.seed.here( fcs.seed.base, "calculate cross-entropy test for tsne" )
  
  ce.diff.test.tsne.res <- ce.diff.test.tsne( 
    dmrd.data, tsne.data, 
    dmrd.event.sample, 
    partition.label = flow.sample.label, 
    partition.color = flow.sample.color, 
    partition.line.type = flow.sample.line.type.single, 
    base.test = fcs.ce.diff.base.test, 
    base.dist = fcs.ce.diff.base.dist, 
    prob.sample.n = fcs.ce.diff.prob.sample.n, 
    dendrogram.order.weight = flow.ce.diff.figure.dendrogram.weight.sample, 
    result = file.path( fcs.ce.diff.tsne.figure.dir, 
                        sprintf( "%s_sample.txt", fcs.ce.diff.tsne.result ) ), 
    cdf.figure = file.path( fcs.ce.diff.tsne.figure.dir, 
                            sprintf( "%s_sample.png", fcs.ce.diff.tsne.figure.cdf ) ), 
    dendrogram.figure = file.path( fcs.ce.diff.tsne.figure.dir, 
                                   sprintf( "%s_sample.png", fcs.ce.diff.tsne.figure.dendrogram ) ) )
  
  
  # calculate cross-entropy test for tsne by cluster---------------
  
  fcs.ce.diff.figure.dendrogram.weight.cluster <- 1 : fcs.cluster.n
  names( fcs.ce.diff.figure.dendrogram.weight.cluster ) <- fcs.cluster
  
  set.seed.here( fcs.seed.base, "calculate cross-entropy test for tsne" )
  
  ce.diff.test.tsne.res <- ce.diff.test.tsne( 
    dmrd.data, tsne.data, 
    dmrd.event.cluster, 
    partition.label = fcs.cluster.label, 
    partition.color = fcs.cluster.color, 
    partition.line.type = fcs.cluster.line.type, 
    base.test = fcs.ce.diff.base.test, 
    base.dist = fcs.ce.diff.base.dist, 
    prob.sample.n = fcs.ce.diff.prob.sample.n, 
    dendrogram.order.weight = fcs.ce.diff.figure.dendrogram.weight.cluster, 
    result = file.path( fcs.ce.diff.tsne.figure.dir, 
                        sprintf( "%s_cluster.txt", fcs.ce.diff.tsne.result ) ), 
    cdf.figure = file.path( fcs.ce.diff.tsne.figure.dir, 
                            sprintf( "%s_cluster.png", fcs.ce.diff.tsne.figure.cdf ) ), 
    dendrogram.figure = file.path( fcs.ce.diff.tsne.figure.dir, 
                                   sprintf( "%s_cluster.png", fcs.ce.diff.tsne.figure.dendrogram ) ) )
  
  # calculate cross-entropy test for umap by condition---------------
  
  set.seed.here( fcs.seed.base, "calculate cross-entropy test for umap" )
  
  ce.diff.test.umap.res <- ce.diff.test.umap( 
    dmrd.knn$dist, dmrd.knn$idx, 
    umap.data, umap.result, 
    dmrd.event.condition, 
    partition.label = fcs.condition.label, 
    partition.color = fcs.condition.color, 
    partition.line.type = fcs.condition.line.type, 
    base.test = fcs.ce.diff.base.test, 
    base.dist = fcs.ce.diff.base.dist, 
    prob.sample.n = fcs.ce.diff.prob.sample.n, 
    dendrogram.order.weight = fcs.ce.diff.figure.dendrogram.weight.condition, 
    result = file.path( fcs.ce.diff.umap.figure.dir, 
                        sprintf( "%s_condition.txt", fcs.ce.diff.umap.result ) ), 
    cdf.figure = file.path( fcs.ce.diff.umap.figure.dir, 
                            sprintf( "%s_condition.png", fcs.ce.diff.umap.figure.cdf ) ), 
    dendrogram.figure = file.path( fcs.ce.diff.umap.figure.dir, 
                                   sprintf( "%s_condition.png", fcs.ce.diff.umap.figure.dendrogram ) ) )
  
  
  # calculate cross-entropy test for umap by sample---------------
  
  set.seed.here( fcs.seed.base, "calculate cross-entropy test for umap" )
  
  ce.diff.test.umap.res <- ce.diff.test.umap( 
    dmrd.knn$dist, dmrd.knn$idx, 
    umap.data, umap.result, 
    dmrd.event.sample, 
    partition.label = flow.sample.label, 
    partition.color = flow.sample.color, 
    partition.line.type = flow.sample.line.type.single, 
    base.test = fcs.ce.diff.base.test, 
    base.dist = fcs.ce.diff.base.dist, 
    prob.sample.n = fcs.ce.diff.prob.sample.n, 
    dendrogram.order.weight = flow.ce.diff.figure.dendrogram.weight.sample, 
    result = file.path( fcs.ce.diff.umap.figure.dir, 
                        sprintf( "%s_sample.txt", fcs.ce.diff.umap.result ) ), 
    cdf.figure = file.path( fcs.ce.diff.umap.figure.dir, 
                            sprintf( "%s_sample.png", fcs.ce.diff.umap.figure.cdf ) ), 
    dendrogram.figure = file.path( fcs.ce.diff.umap.figure.dir, 
                                   sprintf( "%s_sample.png", fcs.ce.diff.umap.figure.dendrogram ) ) )
  
  
  # calculate cross-entropy test for umap by cluster---------------
  
  set.seed.here( fcs.seed.base, "calculate cross-entropy test for umap" )
  
  ce.diff.test.umap.res <- ce.diff.test.umap( 
    dmrd.knn$dist, dmrd.knn$idx, 
    umap.data, umap.result, 
    dmrd.event.cluster, 
    partition.label = fcs.cluster.label, 
    partition.color = fcs.cluster.color, 
    partition.line.type = fcs.cluster.line.type, 
    base.test = fcs.ce.diff.base.test, 
    base.dist = fcs.ce.diff.base.dist, 
    prob.sample.n = fcs.ce.diff.prob.sample.n, 
    dendrogram.order.weight = fcs.ce.diff.figure.dendrogram.weight.cluster, 
    result = file.path( fcs.ce.diff.umap.figure.dir, 
                        sprintf( "%s_cluster.txt", fcs.ce.diff.umap.result ) ), 
    cdf.figure = file.path( fcs.ce.diff.umap.figure.dir, 
                            sprintf( "%s_cluster.png", fcs.ce.diff.umap.figure.cdf ) ), 
    dendrogram.figure = file.path( fcs.ce.diff.umap.figure.dir, 
                                   sprintf( "%s_cluster.png", fcs.ce.diff.umap.figure.dendrogram ) ) )
}



