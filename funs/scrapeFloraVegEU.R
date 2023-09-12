get.LHS.FloraVegEU <- function(species) {
  require(rvest)
  
  url <- paste0('https://floraveg.eu/taxon/overview/',species)
  webpage <- read_html(URLencode(url))
    
  #1. Get plant height ###
  str_plant.height <- '#panel_2 .featureDetail:nth-child(1) b'
  plant.height_html <- html_nodes(webpage, str_plant.height)
  plant.height <- html_text(plant.height_html)
  plant.height <- gsub('[^0-9.-]', '', plant.height)
  if(identical(plant.height, character(0))) {plant.height=NA} else {
    if(plant.height=='' | plant.height=='-') {plant.height=NA} else {
      plant.height <- as.numeric(plant.height)
    }
    }
    
  #2. get specific leaf area###
  str_SLA <- '#panel_3 .align-items-center div'
  SLA_html <- html_nodes(webpage, str_SLA)
  SLA <- html_text(SLA_html)
  if(identical(SLA, character(0))){
    SLA <- NA
   } else {
     for (k in 1:length(SLA)) {
       i2 <- gsub('mm2', '', SLA[[k]])
       SLA[[k]] <- gsub('[^0-9.-]', '', i2)
     }
     SLA <- suppressWarnings(as.numeric(SLA[!is.na(as.numeric(SLA))]))
     
     if (identical(SLA, numeric(0))){
       SLA <- NA
     }
     
  }
    
  #3. get seed mass###
  str_SM <- '#panel_5 b'
  SM_html <- html_nodes(webpage, str_SM)
  SM <- html_text(SM_html)
  SM <- gsub('[^0-9.-]', '', SM)
  SM <- ifelse(identical(SM, character(0)), NA, SM)
  SM <- as.numeric(SM)
    
  res <- data.frame(
      plant_height = plant.height,
      sla = SLA,
      seed_mass = SM)

  return(res)
}

get.LHS.FloraVegEU('Bellis perennis')

#An example of species list (vector):
splis = c('Bellis perennis', #Existing species with full LHS information
          'Pinus mugo', #Existing species with incomplete LHS information
          'Pseudobetckea caucasica', #Existing species with no information of functional traits
          'Taraxacum midolii' #Unexisting species
          )

library(purrr)
get.LHS.FloraVegEU_possibly <- possibly(get.LHS.FloraVegEU, 
                                        data.frame(plant_height = NA, sla = NA, seed_mass = NA))
res <- map(splis, get.LHS.FloraVegEU_possibly) #map the function
names(res) <- splis #set species names
res.df <- do.call(rbind, res) # bind dataset
res.df



#Plot inclomplete cor matrix
# requires tidyverse, reshape2, and corrr packages!!!

#1. Accessory functions ####
getsamplesize <- function(vec_a, vec_b) {
  nrow(drop_na(data.frame(vec_a,vec_b)))
}
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

#2. load the function
plot.inc.cor <- function(mat, round.label.digit=2, cor.method='pearson', size.label=4){
  
  corm <- cor(mat, method = cor.method, use="pairwise.complete.obs")
  
  pcorm <- cor(mat, method = cor.method, use="pairwise.complete.obs") %>%
    get_lower_tri %>%
    reshape2::melt()
  pcorm$value = ifelse(pcorm$Var1 == pcorm$Var2, NA, pcorm$value)
  psams <- corrr::colpair_map(mat, getsamplesize) %>%
    as.data.frame %>%
    tibble::column_to_rownames('term') %>%
    as.matrix %>%
    get_upper_tri %>%
    reshape2::melt() %>%
    setNames(names(pcorm)) 
  
  pcorm$value = round(pcorm$value, round.label.digit)
  # Format the numeric values with two decimal places
  pcorm$formatted_value <- sprintf(paste0('%.',round.label.digit,'f'), pcorm$value)
  pcorm$formatted_value <- ifelse(pcorm$formatted_value=='NA',NA,pcorm$formatted_value)
  
  p <- ggplot() + 
    geom_tile(data = pcorm, aes(x=Var1, y=Var2, fill=value), color='grey80')+
    geom_text(data = psams, aes(x=Var1, y=Var2, label = value, alpha=log(value)), color = "black", size = size.label) +
    geom_text(data = pcorm, aes(x=Var1, y=Var2, label = formatted_value), color = "black", size = size.label) +
    scale_fill_gradient2(low = "#4A6FE3", mid = "white", high = "#D33F6A", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         na.value = 'white',
                         name=paste0(firstup(cor.method),'\ncorrelation')) +
    theme_classic() +
    guides(alpha='none')+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, color='black'),
          axis.text.y = element_text(color='black'),
          axis.title = element_blank(),
          axis.line = element_blank(),
          panel.background = element_rect(fill = "white",
                                          colour = "white")) +
    coord_fixed() +
    geom_abline(slope = 1, intercept = 0, color='grey80') 
  return(p)
}



set.seed(9)

library(tidyverse); library(reshape2); library(corrr)

# Number of rows and columns in the dataframe
num_species <- 100
num_traits <- 7

# Generate random data with NAs
dmatrix <- matrix(runif(num_species * num_traits), nrow = num_species)

# Randomly replace 25% of values with NA
dmatrix[sample(length(dmatrix), length(dmatrix)*0.25)] <- NA

# Convert the matrix to a dataframe
dmatrix <- as.data.frame(dmatrix)
names(dmatrix) = paste0('Trait',1:ncol(dmatrix))

plot.inc.cor(dmatrix, cor.method='pearson', round.label.digit=2, size.label=4)