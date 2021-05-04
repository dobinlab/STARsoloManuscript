grep 'Elapsed (wall clock) time' ../count/salmon-alevin-fry_1.4.0_0.1.0/human_CR_3.0.0/standard/sketch_rad/10X/3/pbmc_5k/*/b*/log.* | awk -f extractTime.awk
