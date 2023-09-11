
get.LHS.FloraVegEU <- function(species) {
  require(rvest)
  cat(species,' |> ')
  url <- paste0('https://floraveg.eu/taxon/overview/',species)
  webpage <- read_html(URLencode(url))
  
  #1. Get plant height ###
  str_plant.height <- '#panel_2 .featureDetail:nth-child(1) b'
  plant.height_html <- html_nodes(webpage, str_plant.height)
  plant.height <- html_text(plant.height_html)
  plant.height <- gsub("[^0-9.-]", "", plant.height)
  if(identical(plant.height, character(0))) {plant.height=NA} else {
    if(plant.height=='' | plant.height=='-') {plant.height=NA}
  }
  
  #2. get specific leaf area###
  str_SLA <- '#panel_3 .align-items-center div'
  SLA_html <- html_nodes(webpage, str_SLA)
  SLA <- html_text(SLA_html)
  if(identical(SLA, character(0))){
    SLA <- NA
  } else {
    if(is_empty(SLA)){
      SLA <- NA
    } else {
      for (i in 1:length(SLA)) {
        i2 <- str_remove_all(SLA[[i]],'mm2')
        SLA[[i]] <- gsub("[^0-9.-]", "", i2)
      }
      SLA <- keep(SLA, function(x){nchar(x)>0})
      if(length(SLA)>1){
        SLA <- keep(SLA, function(x){nchar(x)>1})
      }
      if(SLA=='-'){
        SLA <- NA
      }
    }
  }
  
  #2. get seed mass###
  str_SM <- '#panel_5 b'
  SM_html <- html_nodes(webpage, str_SM)
  SM <- html_text(SM_html)
  SM <- gsub("[^0-9.-]", "", SM)
  SM <- ifelse(identical(SM, character(0)), NA, SM)
  
  HLS = data.frame(plant_height = plant.height,
                   sla = SLA,
                   seed_mass = SM)
  
  return(HLS)
}
