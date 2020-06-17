prelim_plate <- read_excel("Prelim Data/COVID_Methylation_Manifest_v1.xlsx") %>%
  drop_na("Sample ID") %>% mutate(Full_ID = str_c(`Sample ID`, "WG6936921-DNA", sep = "-"))

plate_1_2 <- read_excel("Manifests/COVID_Methylation_Manifest_Plates_1&2_updated.xlsx") %>%
  drop_na("Sample ID") %>% mutate(Full_ID = str_c(`Sample ID`, `Plate Barcode or Number`, sep = "-"))

plate_3_5 <- read_excel("Manifests/COVID_Methylation_Manifest_Plates_3-6_v1_with_BSC_pos.xlsx", sheet = "USED") %>% 
  drop_na("Sample ID") %>% mutate(Full_ID = str_c(`Sample ID`, `Original Plate Barcode`, sep = "-"))

plate_6 <- plate_6 <- read_excel("Manifests/COVID_Methylation_Manifest_Plates_6_replate.xlsx") %>%
  drop_na("Sample ID") %>% mutate(Full_ID = str_c(`Sample ID`, `Original Plate Barcode`, sep = "-"))

prelim_plate %>% select("Covid +/-") %>% mutate_all(as_factor) %>% summary()

plate_1_2 %>% select("Covid +/-") %>% mutate_all(as_factor) %>% summary()

plate_3_5 %>% select("Covid +/-") %>% mutate_all(as_factor) %>% summary()

plate_6 %>% select("Covid +/-") %>% mutate_all(as_factor) %>% summary()

plate7 %>% select("Covid +/-") %>% mutate_all(as_factor) %>% summary()


orig_1_2 <- 