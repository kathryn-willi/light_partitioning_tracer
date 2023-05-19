library(feather)
library(tidyverse)
library(nhdplusTools)
library(sf)
library(mapview)

sf::sf_use_s2(FALSE)

dts <- read_feather('data/simul.feather')

sites <- read_feather('data/simul.feather') %>%
  distinct(SiteID, lat, long) %>%
  # don't know actual CRS but am using NAD83 for comparability to NHD (for now):
  sf::st_as_sf(coords = c("long", "lat"), crs = 4269) %>%
  mutate(comid_tool = NA)

mapview::mapview(sites)
  
# CONTUS flowline features saved data/nhd folder (locally):
path <- list.files("data/nhd", full.names = TRUE)[1]
  
nhd_flowlines <- sf::st_read(dsn = path,
                             layer = "NHDFlowline_Network")

# get the comids for all sites:
# for(i in 1:nrow(sites)){
#   try(sites$comid_tool[i] <- nhdplusTools::discover_nhdplus_id(sites[i,]))
# }
# 
# # where are the NAs?
# nas <- filter(sites, is.na(comid))
# # ... along coasts and lakes. This dataset will require a different workflow... 
# mapview(nas)

# For this dataset it may be better to use `st_nearest()` with CONTUS NHDPlus V2 catchments (downloaded locally):
nhd_catchments <- sf::st_read(dsn = path,
                              layer = "Catchment") %>%
  # add NHD meta data we want to each feature
  left_join(st_drop_geometry(select(nhd_flowlines, comid = COMID, FTYPE, FCODE)), by = c("FEATUREID" = "comid")) %>%
  dplyr::select(comid = FEATUREID,
                FTYPE,
                FCODE) %>%
  # index
  rowid_to_column()

all_comids <- sites %>%
  sf::st_nearest_feature(., nhd_catchments) %>%
  as_tibble()

site_comparison <- cbind(sites, all_comids) %>%
  dplyr::left_join(as_tibble(nhd_catchments), by = c("value" = "rowid")) %>%
  # ID coastal locations that make tracing funky: 
  mutate(COASTAL = ifelse(FTYPE == "Coastline" | FCODE == 56600, "COASTAL", NA)) %>%
  select(SiteID,
         # comid_tool,
         comid,
         COASTAL)

# remove coastal sites from tracing exercise:
sites_noncoastal <- site_comparison %>%
  filter(is.na(COASTAL))

# are there instances where the two methods 
# of getting comids produce different outcomes?
# test <- filter(site_comparison, comid_tool != comid)
# 
# site_catchments <- nhd_catchments %>%
#   dplyr::filter(comid %in% test$comid | comid %in% test$comid_tool)


# ... only where the site's coordinates are sandwiched between
# two comid's. We will use the `st_nearest()` workflow since it can more
# reliably select comid's for all sites
# mapview(test) + mapview(site_catchments)

# back up comid pull (getting all these comid's took a while)
saveRDS(site_comparison, 'data/sites_comid.RDS')

# read in the NHD as a table
nhd <- nhd_flowlines %>% 
  as_tibble() %>%
  dplyr::rename(comid = COMID)

# add nhd meta data to all sites, preserve all data:
site_nhd <- site_comparison %>% 
  as_tibble(.) %>% 
  left_join(nhd, by = "comid")

# function that, for every site, lists other sites up-
# or downstream of it. 
# sites_in_network <- function(network_union){
#   
#   tracer <- function(locations){
#     
#     print(paste0("Tracing ", locations))
#     
#     #site_nhd <- as_tibble(site_nhd)
#     
#     outlet <- site_nhd %>%
#       filter(SiteID == locations)
#     
#     upstream_nhd <- get_UT(nhd, outlet$comid) %>% 
#       as_tibble() %>%
#       rename(comid = value)
#     
#     downstream_nhd <- get_DM(nhd, outlet$comid) %>%
#       as_tibble() %>%
#       rename(comid = value)
#     
#     final <- rbind(upstream_nhd, downstream_nhd) %>%
#       distinct(comid, .keep_all = TRUE) %>%
#       inner_join(., site_comparison, by = 'comid') %>%
#       # All we really want is the SiteIDs up-/down-stream
#       # of the mapped-over site:
#       select(SiteID) 
#     
#     # save each SiteID's trace:
#     write_csv(final, paste0('data/site_trace/', outlet$SiteID, ".csv"))
#     
#     return(final)
#     
#     }
#   
#   crawler <- map_dfr(network_union, tracer)
#   
# }
# 
# networked_sites <- site_comparison[400:2960,] %>%
#   mutate(site_list = map(sites$SiteID, sites_in_network))


tracer <- function(locations){
  
  print(paste0("Tracing ", locations))
  
  outlet <- site_nhd %>%
    filter(SiteID == locations)
  
  upstream_nhd <- get_UT(nhd, outlet$comid, distance = 100 + outlet$LENGTHKM) %>% 
    as_tibble() %>%
    rename(comid = value)
  
  downstream_nhd <- get_DM(nhd, outlet$comid, distance = 100 + outlet$LENGTHKM) %>%
    as_tibble() %>%
    rename(comid = value)
  
  final <- rbind(upstream_nhd, downstream_nhd) %>%
    distinct(comid, .keep_all = TRUE) %>%
    left_join(., site_comparison, by = 'comid') %>%
    mutate(origin = locations)
  
  # save each SiteID's trace:
  write_csv(final, paste0('data/site_trace_100/', str_replace_all(outlet$SiteID, "/","-"), ".csv"))
  
}

sites_noncoastal$SiteID %>%
  walk(~ tracer(locations = .))


# Now that we have lists of every sites' up-/down-stream network and where the related gages are,
# we ccan calculate the distances between each of those sites in each "network":

lengther <- function(sites_related){
  
  file <- read_csv(sites_related) %>%
    #sites_related %>%
    inner_join(select(st_drop_geometry(site_comparison), SiteID, comid), by = c('origin' = 'SiteID')) %>%
    rename(site_1 = origin,
           comid_1 = comid.y,
           site_2 = SiteID,
           comid_2 = comid.x)
    
  print(paste0("Getting combos for ", unique(file$site_1)))
  
  sub_nhd <- nhd_flowlines %>%
    dplyr::filter(COMID %in% file$comid_2)
  
  nhd_prepped <- sub_nhd %>%
    get_tocomid(., add=TRUE) %>%
    mutate(ID=comid, toID=tocomid) #%>%
    dplyr::filter(ID %in% sites$comid_2)
    
  sites <- file %>%
    filter(!is.na(comid_2))
  
  site_nhd <- nhd_prepped %>%
    dplyr::filter(comid %in% sites$comid_2)
  
get_path_lengths(site_nhd$ID, network=nhd_prepped) %>%
  filter(network_distance_km < 100) %>%
    mutate(comid1 = as.numeric(ID_1), comid2 = as.numeric(ID_2)) %>%
    left_join(select(site_comparison, SiteID, comid), by = c("comid1"="comid")) %>%
    rename(site_1=SiteID) %>%
    left_join(select(site_comparison, SiteID, comid), by = c("comid2"="comid")) %>%
    rename(site_2=SiteID) %>%
    select(site_1,site_2,network_distance_km) %>%
    mutate(combos=paste0(site_1,'-',site_2)) %>%
    distinct(.keep_all=TRUE) %>%
    filter(!is.na(site_1),
           !is.na(site_2)) %>%
  
  write_csv(., paste0('data/combo_lengths/', str_replace_all(unique(file$site_1), "/","-"), ".csv"))

  #return(path_lengths)
  
}

# For biggest (43.5 MB) network size, I got this error:
# Error: C stack usage  7953912 is too close to the limit

# sites_related <- read_csv('data/site_trace/MDE_FIELDSERVICES_WQX-XKH4450.csv')
# sites_related <- read_csv('data/site_trace/USGS-07374000.csv')

# list of all sites and the other sites within their flow network:
sites_related <- list.files('data/site_trace/', full.names = TRUE) %>%
  map(~lengther(.)) %>%
  bind_rows() %>%
  mutate(combo1 = paste0(site_1, '-', site_2),
         combo2 = paste0(site_2, '-', site_1))




