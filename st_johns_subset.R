library(feather)
library(tidyverse)
library(nhdplusTools)
library(sf)
library(mapview)

sf::sf_use_s2(FALSE)

dts <- read_feather('data/no_secchi.feather')

sites <- read_feather('data/no_secchi.feather') %>%
  distinct(SiteID, lat, long) %>%
  # don't know actual CRS but am using NAD83 for comparability to NHD (for now):
  sf::st_as_sf(coords = c("long", "lat"), crs = 4269) %>%
  mutate(comid_tool = NA)

mapview::mapview(sites)

# furthest downstream location in St. Johns: 	21FLSJWM-JAXSJR01

# CONTUS flowline features saved data/nhd folder (locally):
path <- list.files("data/nhd", full.names = TRUE)[1]

nhd_flowlines <- sf::st_read(dsn = path,
                             layer = "NHDFlowline_Network")

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

nhd_noncoastal_catchments <- nhd_catchments %>%
  filter(FTYPE != "Coastline" & FCODE != 56600) %>%
  select(-rowid) %>%
  rowid_to_column()

all_comids <- sites %>%
  sf::st_nearest_feature(., nhd_noncoastal_catchments) %>%
  as_tibble()

site_comparison <- cbind(sites, all_comids) %>%
  dplyr::left_join(as_tibble(nhd_noncoastal_catchments), by = c("value" = "rowid")) %>%
  # ID coastal locations that make tracing funky: 
  mutate(COASTAL = ifelse(FTYPE == "Coastline" | FCODE == 56600, "COASTAL", NA)) %>%
  select(SiteID,
         # comid_tool,
         comid,
         COASTAL) #%>%
  # manuall change outlet to be nearest stream instead of a coastal comid:
  # mutate(comid = ifelse(SiteID == "21FLSJWM-JAXSJR01", 16659839, comid))

# remove coastal sites from tracing exercise:
sites_noncoastal <- site_comparison %>%
  filter(is.na(COASTAL))

# read in the NHD as a table
nhd <- nhd_flowlines %>% 
  as_tibble() %>%
  dplyr::rename(comid = COMID)

# add nhd meta data to all sites, preserve all data:
site_nhd <- site_comparison %>% 
  as_tibble(.) %>% 
  left_join(nhd, by = "comid")

tracer <- function(locations){
  
  print(paste0("Tracing ", locations))
  
  outlet <- site_nhd %>%
    filter(SiteID == locations)
  
  upstream_nhd <- get_UT(nhd, outlet$comid) %>% 
    as_tibble() %>%
    rename(comid = value)
  
  downstream_nhd <- get_DM(nhd, outlet$comid) %>%
    as_tibble() %>%
    rename(comid = value)
  
  final <- rbind(upstream_nhd, downstream_nhd) %>%
    distinct(comid, .keep_all = TRUE) %>%
    left_join(., site_comparison, by = 'comid') %>%
    mutate(origin = locations)
  
  return(final)
}

# sites_noncoastal$SiteID %>%
#   walk(~ tracer(locations = .))

# St. Johns's most downstream sample site:
stj <- tracer(locations =  "21FLSJWM-JAXSJR01")
  

# Now that we have lists of St. Johns's up-/down-stream network and where the related sites are,
# we can calculate the distances between each of those sites in each "network":

lengther <- function(sites_related){
  
  file <- #read_csv(sites_related) %>%
    sites_related %>%
    left_join(select(st_drop_geometry(site_comparison), SiteID, comid), by = c('origin' = 'SiteID')) %>%
    rename(site_1 = origin,
           comid_1 = comid.y,
           site_2 = SiteID,
           comid_2 = comid.x)
  
  print(paste0("Getting combos for ", unique(file$site_1)))
  
  sub_nhd <- nhd_flowlines %>%
    dplyr::filter(COMID %in% file$comid_2)
  
  site_nhd <- sub_nhd %>%
    get_tocomid(., add=TRUE) %>%
    mutate(ID=comid, toID=tocomid) #%>%
  #dplyr::filter(ID %in% sites$comid_2)
  
  lengths <- get_path_lengths(site_nhd$ID, network=nhd_prepped) %>%
    # make sure we still capture sites on the same comid as furthest-down 
    # point:
    rbind(tibble(ID_1 = as.character(filter(nhd_prepped, toID == 0)$comid),
                 ID_2 = as.character(filter(nhd_prepped, toID == 0)$comid),
                 network_distance_km = 0)) %>%
    #filter(network_distance_km < 100) %>%
    mutate(comid1 = as.numeric(ID_1), comid2 = as.numeric(ID_2)) %>%
    inner_join(select(site_comparison, SiteID, comid), by = c("comid1"="comid")) %>%
    rename(site_1=SiteID) %>%
    inner_join(select(site_comparison, SiteID, comid), by = c("comid2"="comid")) %>%
    rename(site_2=SiteID) %>%
    select(site_1,site_2,network_distance_km) %>%
    mutate(combos=paste0(site_1,'-',site_2)) %>%
    distinct(.keep_all=TRUE) %>%
    filter(!is.na(site_1),
           !is.na(site_2)) %>%
    #filter(grepl(unique(file$site_1), site_1)) %>%
    distinct(combos,.keep_all = TRUE)
    
   # write_csv(., paste0('data/combo_lengths/', str_replace_all(unique(file$site_1), "/","-"), ".csv"))
  
  return(lengths)
  
}

st_j_distances <- lengther(sites_related = stj) 

st_j_mapper <- st_j_distances %>%
  filter(grepl(unique(file$site_1), combos))

full_org <- st_j_mapper %>%
  filter(site_1 != unique(file$site_1)) %>%
  select(SiteID = site_1,
         network_distance_km)

full_org_2 <- st_j_mapper %>%
  filter(site_2 != unique(file$site_1)) %>%
  select(SiteID = site_2,
         network_distance_km)

fullest_org <- rbind(full_org, full_org_2) %>%
  distinct(.keep_all = TRUE)

for_matt <- site_comparison %>%
  right_join(.,fullest_org, by="SiteID") %>%
  select(-COASTAL) %>%
  mutate(origin = unique(file$site_1))
 
mapview(for_matt, zcol = "network_distance_km") + mapview(filter(nhd_flowlines, COMID %in% site_nhd$ID)) +mapview(site_comparison)
