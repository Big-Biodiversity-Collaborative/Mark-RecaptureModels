# Plot RMNP banding sites and BTLH wintering sightings in Mexico
# Gaby Samaniego
# gabysamaniego@arizona.edu
# 2025-03-23

# Load packages
library(tidyverse)
library(sf)
library(ggrepel)
library(patchwork)

# Clear environment
rm(list = ls())

# ------------------------------- LOAD DATA ---------------------------------- #

# Load BTLH range map (source eBird) 
range <- st_read('data/sites-BTLH-range-map/brthum_range_2021.gpkg')

# Subset breeding and non breeding ranges 
breeding <- range %>% 
  filter(season == 'breeding')

nonbreeding <- range %>% 
  filter(season == 'nonbreeding')

# Read RMNP banding sites and exclude sites that won't be part of analysis  
summer.sites <- read.csv('data/sites-BTLH-range-map/thinned-summer-sites.csv') %>% 
  filter(!site %in% c('WB2','WB1', 'WPK1', 'NFPC', 'POLC', 'SHIP'))

# Convert data frame to sf object and add column useful to plot the data
summer.sf <- st_as_sf(summer.sites, 
                      coords = c('longitude', 'latitude'), 
                      crs = st_crs(breeding)) %>%
  mutate(point_type = 'BTLH banding sites')

# Read Mexico's sites, convert to sf and add column useful to plot the data
winter.sites <- read.csv('data/sites-BTLH-range-map/thinned-winter-sites.csv')
winter.sf <- st_as_sf(winter.sites, 
                      coords = c('longitude', 'latitude'), 
                      crs = st_crs(nonbreeding)) %>%
  mutate(point_type = 'BTLH sightings')

# RMNP boundary
rmnp.sf <- st_read('data/sites-BTLH-range-map/RMNP-boundary/Boundary_(Line).shp')

# ---------------------- PLOT SUMMER AND WINTERING GROUNDS ------------------- # 
# ---------------- INCLUDING RMNP BOUNDARY AND MEXICO SIGHTINGS -------------- #

# Create dummy data to force fill legend in the map
legend.fill <- data.frame(x = c(0, 0),
                          y = c(0, 0),
                          range_type = c('Breeding', 'Nonbreeding'))

# Create map
main.map <- ggplot() + 
  
  # Plot breeding and non breeding range 
  geom_sf(data = breeding, 
          fill = '#f6b48f', # Pale orange, color blind friendly
          color = NA) +
  geom_sf(data = nonbreeding, 
          fill = '#a6d9ce', # Pale teal, color blind friendly
          color = NA) +
  
  # Add RMNP boundary
  geom_sf(data = rmnp.sf,
          aes(color = 'Rocky Mountain National Park Boundary'),
          fill = NA,
          size = 0.5) +
  
  # Plot borders
  borders('world', colour = 'gray18') +
  borders('state', colour = 'gray18') +
  
  # Plot GBIF sightings
  geom_sf(data = winter.sf, 
          aes(shape = point_type), 
          size = 1.5, 
          color = 'gray18') +
  
  # Create invisible dummy tiles to create fill legend
  geom_tile(data = legend.fill, 
            aes(x = x, y = y, fill = range_type),
            alpha = 0) +
  
  # Manually fill colors for breeding and non breeding
  scale_fill_manual(values = c('Breeding' = '#f6b48f', 
                               'Nonbreeding' = '#a6d9ce'),
                    name = 'Broad-tailed Hummingbird Range',
                    labels = c('Breeding' = 'Summer Grounds', 
                               'Nonbreeding' = 'Wintering Grounds')) +
  
  # Add color scale for RMNP outline
  scale_color_manual(values = c('Rocky Mountain National Park Boundary' = 'gray18'),
                     name = '') +
  
  # Define point's shapes
  scale_shape_manual(values = c('BTLH sightings' = 17),     # Triangle
                     name = '',
                     labels = c('BTLH sightings' = 'GBIF sightings')) +
  
  # Define map extent
  coord_sf(expand = FALSE, 
           xlim = c(-118, -89),  
           ylim = c(13, 42.5)) +
  
  # Edit map's theme
  theme_bw(base_size = 7) +
  theme(panel.grid.major = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x  = element_text(size = 7),
        axis.text.y  = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.text  = element_text(size = 7),
        legend.key.size = unit(0.2, 'cm'),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.box.background = element_blank(),
        legend.box = 'vertical',
        legend.spacing.y = unit(-3, 'mm'),
        legend.position = c(0.29, 0.12)) +  
  
  # Customize legend appearance and order
  guides(fill = guide_legend(order = 1,
                             override.aes = list(size = 3, alpha = 1)),
         shape = guide_legend(order = 2,
                              override.aes = list(size = 2)),
         color = guide_legend(order = 3))
main.map

# Save plot
ggsave(path = 'output/plots',
       filename = 'BTHU-range-GBIF-sightings.png',
       plot = main.map,
       device = 'png',
       dpi = 300,
       width = 6,
       height = 6.5)

# --------------------------- PLOT RMNP BANDING SITES ------------------------ # 
# -------------------------- INCLUDING PARK'S BOUNDARY ----------------------- #

# Plot banding sites inside RMNP
RMNP.map <- ggplot() + 
  
  # Add banding sites
  geom_sf(data = summer.sf,
          aes(shape = 'Banding Site'),
          size = 1.5,
          color = 'gray18') +
  
  # Add RMNP boundary
  geom_sf(data = rmnp.sf,
          aes(color = 'Rocky Mountain National Park Boundary'),
          fill = NA,
          size = 0.5) +
  
  # Set RMNP boundary color to lighter grey
  scale_color_manual(name = '',
                     values = c('Rocky Mountain National Park Boundary' = 'grey60')) +
  
  # Add site labels
  geom_text_repel(data = summer.sf,
                  aes(label = site, geometry = geometry),
                  stat = 'sf_coordinates',
                  size = 2,
                  color = 'gray18',  
                  
                  # Tighter label placement
                  box.padding   = unit(0.6, 'lines'),
                  point.padding = unit(0.5, 'lines'),
                  max.overlaps  = Inf,
                  
                  # Connector lines
                  segment.color = 'gray18', 
                  segment.size  = 0.3,
                  min.segment.length = 0,
                  segment.curvature = 0) +
  
  # Map extent
  coord_sf(expand = FALSE,
           xlim = c(-106, -105.4),
           ylim = c(40.1, 40.6)) +
  
  # Shape scale
  scale_shape_manual(name = '',
                     values = c('Banding Site' = 16)) +
  
  # Fix x axis breaks
  scale_x_continuous(breaks = c(-105.9, -105.8, -105.7, -105.6, -105.5),
                     labels = c('105.9°W', '105.8°W', '105.7°W', '105.6°W', '105.5°W')) +
  
  # Theme
  theme_bw(base_size = 7) +
  theme(panel.grid.major = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x  = element_text(size = 7),
        axis.text.y  = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.text  = element_text(size = 7),
        legend.key.size = unit(0.2, 'cm'),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.box.background = element_blank(),
        legend.box = 'vertical',
        legend.spacing.y = unit(-3, 'mm'),
        legend.position = c(0.35, 0.09))
RMNP.map

# Save plot
ggsave(path = 'output/plots',
       filename = 'RMNP-banding-sites.png',
       plot = RMNP.map,
       device = 'png',
       dpi = 300,
       width = 8,
       height = 5)

# ------------------- COMBINE MAPS IN A TWO PANEL FIGURE --------------------- #

# Combine maps
combined.map <- main.map + RMNP.map +
  plot_layout(ncol = 2, widths = c(1.6, 1.1)) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 10),

    # spacing between panel A and B
    panel.spacing = unit(-0.3, 'cm'),
    
    # outer margin of the full figure
    plot.margin = margin(2, 2, 2, 1))

combined.map

# Save final figure
ggsave(path = 'output/plots',
       filename = 'study-sites-two-panel.png',
       plot = combined.map,
       device = 'png',
       dpi = 400,
       width = 8,   
       height = 4.5)
