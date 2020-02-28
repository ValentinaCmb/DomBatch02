# Packages ----

library("addinslist")
library(tidyverse)
library(readr)
library(RColorBrewer)
library(modelr)
library(kableExtra)
library(ggThemeAssist)


getwd()

# LS = Lipid Search Software i.e. table of lipids etc obtained from LipidSearch
# CD =  Compund Discoverer

##%######################################################%##
#                                                          #
####                   1st function                     ####
#                                                          #
##%######################################################%##


#  load and transorm the data obtained in LipidSearch in a format can work on in R 

readLipidSearch <- function(filename) {
  skiprow <- read_tsv(filename)
  skiprow <- which(grepl("Result", skiprow$`Job Name`)) +3
  df <- read_tsv(filename, skip = skiprow) %>%
    separate(`FA Group Key`, c("TC", "DB"), sep = ":", remove = FALSE)
  }


LSdata <- readLipidSearch("./Data/LS_TableDomBatch01.txt")
#saved in Z:\Documents\Lipidomics Domestication Batch01\R\Data


ggplot(LSdata, aes(BaseRt, `Calc Mass`, colour = Class)) +geom_point() 

# first check to see if some compounds look a bit "dodgy"; so separate the different lipid classes
# but want to look to see if there are any difference

ggplot(filter(LSdata, Class == "TG"), aes(BaseRt, `Calc Mass`, colour = DB)) +
  geom_point() +
  facet_wrap(~Class)


##%######################################################%##
#                                                          #
####                    2nd function                    ####
#                                                          #
##%######################################################%##
 

# following step 1, now we want to pull out  data from the lipidsearch table to make a new table that will run in CD

LS_to_CD_list <- function(df){
  #unassigned <- filter(df, Rej == 0) #need to change to do this only if there is a column called Rej
  calc <- df %>%
    ungroup() %>%
    select(Name = LipidMolec, LSFormula = Formula, RT = BaseRt, MainIon, `Calc Mass`, Class, FA, `FA Group Key`) %>%
    unique() %>%
    mutate(C = as.numeric(str_match(LSFormula, "\\bC(\\d+)")[,2]),
           H = as.numeric(str_match(LSFormula, "\\bH(\\d+)")[,2]),
           N = as.numeric(ifelse(grepl("N", LSFormula), 
                                 ifelse(grepl("N\\d+", LSFormula), 
                                        str_match(LSFormula, "\\bN(\\d+)")[,2], 1), 0)), 
           O = as.numeric(ifelse(grepl("O", LSFormula), 
                                 ifelse(grepl("O\\d+", LSFormula), 
                                        str_match(LSFormula, "\\bO(\\d+)")[,2], 1), 0)),
           P = as.numeric(ifelse(grepl("P", LSFormula), 
                                 ifelse(grepl("P\\d+", LSFormula), 
                                        str_match(LSFormula, "\\bP(\\d+)")[,2], 1), 0)),
           S = ifelse(grepl("S", LSFormula), 
                      ifelse(grepl("S\\d+", LSFormula), 
                             str_match(LSFormula, "\\bS(\\d+)")[,2], 1), ""),
           N = ifelse(MainIon == "+NH4", N + 1, N),
           H = ifelse(MainIon == "+NH4", H + 3, H),
           Name = ifelse(MainIon == "+NH4", paste(Name, "+NH4", sep = " "), Name),
           Formula = paste("C", C, "H", H, ifelse(N > 0, "N", ""), ifelse(N > 1, N, ""), 
                           "O", O, ifelse(P > 0, "P", ""), ifelse(P > 1, P, ""), sep = "")) 
  # %>%
    # mutate(Lipid_Class = ifelse(grepl("NO6", Formula) & !grepl("P", Formula), "TG", 
                               # ifelse(grepl("N2O5P", Formula), "CAEP", NA)),
           # Lipid_Class = ifelse(grepl("NO8P", Formula) & C %% 2 == 0, "PC" ,
               #                 ifelse(grepl("NO8P", Formula) & C %% 2 != 0, "PE", Lipid_Class)),
           # Lipid_Class = ifelse(C>30 & grepl("NO7P", Formula) & C %% 2 == 0, "PC",
                              #  ifelse(grepl("NO7P", Formula) & C %% 2 != 0, "PE", Lipid_Class)),
          # Lipid_Species = ifelse(grepl("NO6", Formula) & !grepl("P", Formula),
                                 # paste(C-3, ((C-3)*2-(H-5))/2, sep = ":"), NA)) 
  
  
  return(calc)
}

LSdata01 <- LS_to_CD_list(LSdata)
write_csv(LSdata01, path = "Data/Table_CD28Nov.csv")

# this below is not needed 
T01 <- LSdata01 %>% 
      select(Name, LSFormula, RT, MainIon, `Calc Mass`, Class, FA, `FA Group Key`) %>% 
      filter(Class %in% c("DG", "FA", "PE", "LPC", "LPE", "LPI", "LPS", "PA", "PC", "PG", "PI", "PS", "TG"))



##%######################################################%##
#                                                          #
####                      3rd part                      ####
#                                                          #
##%######################################################%##


# now that CD has run we want to plot this data and see what makes sense


CDdata <- read_tsv("Data/CompDiscoBatch01Results11Dec.csv", col_types = cols(Name = col_character()))

key <- read_csv("Data/Table_CD28Nov.csv", col_types = cols(Name = col_character())) %>%  
  select(Name, Class, FA, `FA Group Key`) %>% 
  unique()
  # rename(LSFormula= Formula) %>% 
  # mutate(`FA Group Key` = str_remove(str_replace(`FA Group Key`, ":", "_"), ":00"))
  # 

CDdata_gathered06 <- CDdata %>% 
  mutate(Name = gsub("Oleic acid", "FA(18:1)", x = Name)) %>% ###use quotes not the other thing (back thingy``)
  gather(Sample, Area, contains("Area: ")) %>%
  group_by(Name, Formula, `Molecular Weight`, `RT [min]`) %>%
  filter(!is.na(Area),
         `RT [min]` < 25,
         !is.na(Name)) %>% 
  mutate(count = n()) %>% 
  left_join(key) %>% 
  mutate(Bond = str_extract(`FA Group Key`, "[a-zA-Z]+")) %>% 
  mutate(TC = word(`FA Group Key`, 1, 1, sep = ":"),
         DB = word(`FA Group Key`, -1, -1, sep = ":")) %>% 
  mutate(DB = as.numeric(str_extract(DB, "\\d+"))) %>% 
  group_by(Sample) %>%
  mutate(Area2 = Area/sum(Area) * 1000000) %>% 
  ungroup() %>% 
  group_by(Formula, `Molecular Weight`, `RT [min]`) %>%
  mutate(reference = median(Area2),
         quotient = Area2/reference) %>%
  ungroup() %>%
  group_by(Sample) %>%
  mutate(quotient.median = median(quotient),
         pqn = Area2/quotient.median) # %>%
  # mutate (pop = (str_extract(Sample, "\\w+_\\w+"))) %>% 
  # separate(pop, c("population", "sample"), sep="(?<=[A-Za-z])(?=[0-9])")
  # 
 
CDdata_gathered06 
write.csv(CDdata_gathered06, "./CDdata_gathered06.csv", row.names =  FALSE)


lipids.to.keep <- c("DG", "TG", "PC", "PI", "PA", "PG", "PS", "SM", "LPC", "FA","PE", "Cer")

res <-  CDdata_gathered06 %>% 
  filter(grepl(pattern = paste(lipids.to.keep, collapse = "|"), Name)) %>% 
  filter(!grepl("Similar", Name))
           
Batch01 <- res %>% 
  select(Name, Formula, Class, `FA Group Key`, TC, DB, count, `Molecular Weight`, Sample, `RT [min]`, Area, pqn, Area2)

write_csv(Batch01, path = "Data/Batch01.csv")
view(Batch01)

# ____ THESE 2 BELOW are NOT necessarily NEEDED, only if you want to focus on certauinn classes, otherwise skip _____#
PC <- CDdata_gathered %>% 
  filter(Class == "PC")

ggplot(PC, aes(`RT [min]`, `Molecular Weight`, colour = DB)) +
  geom_point()
#______________________________________________________________ #


# Pre-final adjustment (still need to remove NA)

p1 <- Batch01 %>% 
  group_by(Class) %>% 
  rename( `Lipid species` = `FA Group Key`) %>% 
  filter(!grepl("e|p|t|Q|d|O", `Lipid species`))

write_csv(p1, path = "Data/p1.csv")

# PLOTS: 
# First plot with areas, not conc yet ##
# this gives the boxplot w title "Batch01" 

batch01.plot <-  ggplot(p1, aes(`Lipid species`, log(Area))) +
  geom_boxplot(aes(fill= `Lipid species`)) +
  facet_wrap(~Class, scales ="free")

batch01.plot +
  theme(legend.position ="none", 
        axis.text.x = element_text(angle = 45, vjust =0.5)) +
  labs( title = paste("Batch01"))

#### plot normalised data 
batch01.plot.pqn <-  ggplot(p1, aes(`Lipid species`, log(pqn))) +
  geom_boxplot(aes(fill= `Lipid species`)) +
  facet_wrap(~Class, scales ="free")

batch01.plot.pqn +
  theme(legend.position ="none", 
        axis.text.x = element_text(angle = 45, vjust =0.5)) +
  labs( title = paste("Batch01 - normalised data"))
                                                          
##%######################################################%##
#                                                          #
#    Join the table (Batch01StandardsNew) containing      #
#       the STANDARDS for Batch01                          #
#        (ran with the compounds in the LCMS)              #  
#                 with  Batch01Results                     #
#          i.e. to the compunds of Batch01                 #
#                                                          #
##%######################################################%##

Batch01Standards <- read_csv("./Data/Batch01StandardsNew.csv")
view(Batch01Standards)

batch01.standards <- Batch01Standards %>% 
  select(Class, StandardsName, Conc, slope, intercept)


# calculate the concentrations of the compounds in Batch01 based on the standards in Batch01

batch01.values <-  p1 %>% 
  left_join(batch01.standards, by = "Class")  %>% 
  mutate(conc_mgmL_compounds = as.numeric(exp((log(Area) - intercept)/ slope))) %>% 
  select(-slope, - intercept) %>% 
  filter(!Class %in% c("LPE", "LPI", "LPS", "PA")) %>% 
  filter(!is.na(Class))

#check elimination of those lipids with unique
unique(batch01.values$Class)

batch01.values


# split the Sample column in population w domestication and sample number

batch01.pop <- batch01.values %>% 
        mutate (pop = (str_extract(Sample, "\\w+_\\w+"))) %>% 
        separate(pop, c("population", "sample"), sep="(?<=[A-Za-z])(?=[0-9])")
  
 
unique(batch01.pop$Class)

#str_extract("AB:C_N01.DE", "\\w+_\\w+" )
#unlist(str_split("pop","(?<=[A-Za-z])(?=[0-9])"))
# str_extract("AB:C_N01.DE", "\\d+" )
# str_extract("C_N01", "\\d+" )
# str_remove_all("C_N01", "[[:punct:]]" )
# str_split(str_split("C_N01","_")[[1]][[2]], "(?<=[A-Za-z])(?=[0-9])")
# unlist(str_split("C_N01","(?<=[A-Za-z])(?=[0-9])"))[2]


##%######################################################%##
#                                                          #
#               PLOTS of the CONCENTRATIONS             ####
#                                                         #
##%######################################################%##

# ! batch01.pop has the population now


batch01.plot.conc <-  ggplot(batch01.pop, aes(`Lipid species`, log(conc_mgmL_compounds))) + 
  geom_boxplot(aes(fill= `population`)) +
  facet_wrap(~Class, scales ="free") + 
  theme(plot.subtitle = element_text(vjust = 1), 
    plot.caption = element_text(vjust = 1), 
    axis.text = element_text(size = 7, vjust = -0.25), 
    axis.text.x = element_text(angle = 40), 
    panel.background = element_rect(fill = "palegoldenrod"), 
    legend.position = "right") +labs(title = "Batch01: lipid species concentrations", 
    y = "log concentrations") 

batch01.plot.conc



###normalised area results -  NOT for concentration

batch01.plot.conc.pqn <-  ggplot(batch01.pop, aes(`Lipid species`, log(pqn))) + 
  geom_boxplot(aes(fill= `population`)) +
  facet_wrap(~Class, scales ="free") + 
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text = element_text(size = 7, vjust = -0.25), 
        axis.text.x = element_text(angle = 40), 
        panel.background = element_rect(fill = "palegoldenrod"), 
        legend.position = "right") +labs(title = "Batch01: lipid species normalised areas", 
                                         y = "log(pqn)") 

batch01.plot.conc.pqn
write_csv(batch01.plot.conc.pqn, path = "Data/batch01.plot.conc.pqn.csv")


FAnorm.area <- batch01.pop %>% 
  filter(Class == "FA")

FA.batch01.plot.conc.pqn <-  ggplot(FAnorm.area, aes(`Lipid species`, log(pqn))) + 
  geom_boxplot(aes(fill= `population`)) +
  facet_wrap(~Class, scales ="free") + 
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text = element_text(size = 7, vjust = -0.25), 
        axis.text.x = element_text(angle = 40), 
        #panel.background = element_rect(fill = "palegoldenrod"), 
        legend.position = "right") +labs(title = "Batch01: FA lipid species normalised areas", 
                                         y = "log(pqn)") 

FA.batch01.plot.conc.pqn


PGnorm.area <- batch01.pop %>% 
  filter(Class == "PG")

PG.batch01.plot.conc.pqn <-  ggplot(PGnorm.area, aes(`Lipid species`, log(pqn))) + 
  geom_boxplot(aes(fill= `population`)) +
  facet_wrap(~Class, scales ="free") + 
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text = element_text(size = 7, vjust = -0.25), 
        axis.text.x = element_text(angle = 40), 
        #panel.background = element_rect(fill = "palegoldenrod"), 
        legend.position = "right") +labs(title = "Batch01: PG lipid species normalised areas", 
                                         y = "log(pqn)") 

PG.batch01.plot.conc.pqn


TGnorm.area <- batch01.pop %>% 
  filter(Class == "TG")

TG.batch01.plot.conc.pqn <-  ggplot(TGnorm.area, aes(`Lipid species`, log(pqn))) + 
  geom_boxplot(aes(fill= `population`)) +
  facet_wrap(~Class, scales ="free") + 
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text = element_text(size = 7, vjust = -0.25), 
        axis.text.x = element_text(angle = 40), 
        #panel.background = element_rect(fill = "palegoldenrod"), 
        legend.position = "right") +labs(title = "Batch01: TG lipid species normalised areas", 
                                         y = "log(pqn)") 

TG.batch01.plot.conc.pqn


DGnorm.area <- batch01.pop %>% 
  filter(Class == "DG")

DG.batch01.plot.conc.pqn <-  ggplot(DGnorm.area, aes(`Lipid species`, log(pqn))) + 
  geom_boxplot(aes(fill= `population`)) +
  facet_wrap(~Class, scales ="free") + 
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text = element_text(size = 7, vjust = -0.25), 
        axis.text.x = element_text(angle = 40), 
        #panel.background = element_rect(fill = "palegoldenrod"), 
        legend.position = "right") +labs(title = "Batch01: DG lipid species normalised areas", 
                                         y = "log(pqn)") 

DG.batch01.plot.conc.pqn


LPCnorm.area <- batch01.pop %>% 
  filter(Class == "LPC")

LPC.batch01.plot.conc.pqn <-  ggplot(LPCnorm.area, aes(`Lipid species`, log(pqn))) + 
  geom_boxplot(aes(fill= `population`)) +
  facet_wrap(~Class, scales ="free") + 
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text = element_text(size = 7, vjust = -0.25), 
        axis.text.x = element_text(angle = 40), 
        #panel.background = element_rect(fill = "palegoldenrod"), 
        legend.position = "right") +labs(title = "Batch01: LPC lipid species normalised areas", 
                                         y = "log(pqn)") 

LPC.batch01.plot.conc.pqn


PCnorm.area <- batch01.pop %>% 
  filter(Class == "PC")

PC.batch01.plot.conc.pqn <-  ggplot(PCnorm.area, aes(`Lipid species`, log(pqn))) + 
  geom_boxplot(aes(fill= `population`)) +
  facet_wrap(~Class, scales ="free") + 
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text = element_text(size = 7, vjust = -0.25), 
        axis.text.x = element_text(angle = 40), 
        #panel.background = element_rect(fill = "palegoldenrod"), 
        legend.position = "right") +labs(title = "Batch01: PC lipid species normalised areas", 
                                         y = "log(pqn)") 

PC.batch01.plot.conc.pqn


PEnorm.area <- batch01.pop %>% 
  filter(Class == "PE")

PE.batch01.plot.conc.pqn <-  ggplot(PEnorm.area, aes(`Lipid species`, log(pqn))) + 
  geom_boxplot(aes(fill= `population`)) +
  facet_wrap(~Class, scales ="free") + 
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text = element_text(size = 7, vjust = -0.25), 
        axis.text.x = element_text(angle = 40), 
        #panel.background = element_rect(fill = "palegoldenrod"), 
        legend.position = "right") +labs(title = "Batch01: PE lipid species normalised areas", 
                                         y = "log(pqn)") 

PE.batch01.plot.conc.pqn


PInorm.area <- batch01.pop %>% 
  filter(Class == "PI")

PI.batch01.plot.conc.pqn <-  ggplot(PInorm.area, aes(`Lipid species`, log(pqn))) + 
  geom_boxplot(aes(fill= `population`)) +
  facet_wrap(~Class, scales ="free") + 
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text = element_text(size = 7, vjust = -0.25), 
        axis.text.x = element_text(angle = 40), 
        #panel.background = element_rect(fill = "palegoldenrod"), 
        legend.position = "right") +labs(title = "Batch01: PI lipid species normalised areas", 
                                         y = "log(pqn)") 

PI.batch01.plot.conc.pqn


PSnorm.area <- batch01.pop %>% 
  filter(Class == "PS")

PS.batch01.plot.conc.pqn <-  ggplot(PSnorm.area, aes(`Lipid species`, log(pqn))) + 
  geom_boxplot(aes(fill= `population`)) +
  facet_wrap(~Class, scales ="free") + 
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text = element_text(size = 7, vjust = -0.25), 
        axis.text.x = element_text(angle = 40), 
        #panel.background = element_rect(fill = "palegoldenrod"), 
        legend.position = "right") +labs(title = "Batch01: PS lipid species normalised areas", 
                                         y = "log(pqn)") 

PS.batch01.plot.conc.pqn




############################################################
############################################################
##%######################################################%##
##                                                        ##
##            *  SINGLE Plots  *                          ##
##                                                        ##
##%######################################################%##
############################################################


unique(batch01.pop$Class)
#"PE"  "LPC" "DG"  "PS"  "PI"  "PG"  "TG"  "PC" 


##%######################################################%##
#                                                          #
####                         FA                         ####
#                                                          #
##%######################################################%##

batch01.FA <-  batch01.pop %>% 
  filter(Class == "FA") 

batch01.plot.FA <-  ggplot(batch01.FA, aes(`Lipid species`, log(conc_mgmL_compounds), fill = `population`)) +
  geom_boxplot(width =1) +
  scale_fill_brewer(palette = "Set1") +
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text = element_text(size = 7, vjust = -0.25), 
        axis.text.x = element_text(angle = 40), 
        panel.background = element_rect(fill = "grey91"), 
        legend.position = "right") +
  labs(title = "Batch01 FA - lipid species concentrations per population", y = "Log concentrations", x = "Lipid species (Number of carbon _ Number of double bonds)") + 
  theme(axis.text = element_text(size = 8, face = "bold"), 
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10), 
        legend.text = element_text(size = 11), 
        legend.title = element_text(size = 11))  


batch01.plot.FA
ggsave("batch01.plot.FA.png", width =10, height = 8)
###normalised FFA
batch01.plot.FA.pqn <-  ggplot(batch01.FA, aes(`Lipid species`, log(pqn), fill = `population`)) +
  geom_boxplot(width =1) +
  scale_fill_brewer(palette = "Set1") +
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text = element_text(size = 7, vjust = -0.25), 
        axis.text.x = element_text(angle = 40), 
        panel.background = element_rect(fill = "grey91"), 
        legend.position = "right") +
  labs(title = "Batch01 FA - lipid species normalised area per population", y = "log(pqn)", x = "Lipid species (Number of carbon _ Number of double bonds)") + 
  theme(axis.text = element_text(size = 8, face = "bold"), 
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10), 
        legend.text = element_text(size = 11), 
        legend.title = element_text(size = 11))
batch01.plot.FA.pqn

##%######################################################%##
#                                                          #
####                         PE                         ####
#                                                          #
##%######################################################%##

batch01.PE <-  batch01.pop %>% 
  filter(Class == "PE") 

batch01.plot.PE <-  ggplot(batch01.PE, aes(paste(`Lipid species`, Name), log(conc_mgmL_compounds), fill = `population`)) +
  geom_boxplot(width =1) +
  scale_fill_brewer(palette = "Set1") +
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text = element_text(size = 7, vjust = -0.25), 
        axis.text.x = element_text(angle = 40), 
        panel.background = element_rect(fill = "grey91"), 
        legend.position = "right") +
  labs(title = "Batch01 PE - lipid species concentrations per population", y = "Log concentrations", x = "Lipid species (Number of carbon _ Number of double bonds)") + 
  theme(axis.text = element_text(size = 8, face = "bold"), 
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10), 
        legend.text = element_text(size = 11), 
        legend.title = element_text(size = 11))  


batch01.plot.PE
##normalised area plot
batch01.plot.PE.pqn <-  ggplot(batch01.PE, aes(paste(`Lipid species`, Name), log(pqn), fill = `population`)) +
  geom_boxplot(width =1) +
  scale_fill_brewer(palette = "Set1") +
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text = element_text(size = 7, vjust = -0.25), 
        axis.text.x = element_text(angle = 40), 
        panel.background = element_rect(fill = "grey91"), 
        legend.position = "right") +
  labs(title = "Batch01 PE - lipid species normalised area per population", y = "log(pqn)", x = "Lipid species (Number of carbon _ Number of double bonds)") + 
  theme(axis.text = element_text(size = 8, face = "bold"), 
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10), 
        legend.text = element_text(size = 11), 
        legend.title = element_text(size = 11))  


batch01.plot.PE.pqn

ggsave("batch01.plot.PE.png", width =10, height = 8)

##%######################################################%##
#                                                          #
####                         LPC                        ####
#                                                          #
##%######################################################%##

batch01.pop$Class

batch01.LPC <-  batch01.pop %>% 
  filter(Class == "LPC") 

batch01.plot.LPC <-  ggplot(batch01.LPC, aes(`Lipid species`, log(conc_mgmL_compounds), fill = `population`)) +
  geom_boxplot(width =1) +
  scale_fill_brewer(palette = "Set1") +
    theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text = element_text(size = 7, vjust = -0.25), 
        axis.text.x = element_text(angle = 40), 
        panel.background = element_rect(fill = "grey91"), 
        legend.position = "right") +
        labs(title = "Batch01 LPC - lipid species concentrations per population", y = "Log concentrations", x = "Lipid species (Number of carbon _ Number of double bonds)") + 
    theme(axis.text = element_text(size = 8, face = "bold"), 
    axis.text.x = element_text(size = 10), 
    axis.text.y = element_text(size = 10), 
    legend.text = element_text(size = 11), 
    legend.title = element_text(size = 11))  
    

batch01.plot.LPC
ggsave("batch01.plot.LPC.png", width =10, height = 8)


batch01.LPC <-  batch01.pop %>% 
  filter(Class == "LPC") 

batch01.plot.LPC.pqn <-  ggplot(batch01.LPC, aes(`Lipid species`, log(pqn), fill = `population`)) +
  geom_boxplot(width =1) +
  scale_fill_brewer(palette = "Set1") +
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text = element_text(size = 7, vjust = -0.25), 
        axis.text.x = element_text(angle = 40), 
        panel.background = element_rect(fill = "grey91"), 
        legend.position = "right") +
  labs(title = "Batch01 LPC - lipid species normalised area per population", y = "Log(pqn)", x = "Lipid species (Number of carbon _ Number of double bonds)") + 
  theme(axis.text = element_text(size = 8, face = "bold"), 
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10), 
        legend.text = element_text(size = 11), 
        legend.title = element_text(size = 11))  


batch01.plot.LPC.pqn

##%######################################################%##
#                                                          #
####                      DG plot                       ####
#                                                          #
##%######################################################%##

batch01.DG <-  batch01.pop %>% 
  filter(Class == "DG") %>% 
  mutate(`Lipid species` = gsub("_0", ":", x = `Lipid species`)) ##add this earlier in your script!!

batch01.plot.DG <-  ggplot(batch01.pop %>% filter(Class == "DG"), aes(paste(`Lipid species`, Name, sep = " "), log(conc_mgmL_compounds), fill = `population`)) +
  geom_boxplot(width =1) +
  scale_fill_brewer(palette = "Set1") +
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text = element_text(size = 7, vjust = -0.25), 
        axis.text.x = element_text(angle = 40), 
        panel.background = element_rect(fill = "grey91"), 
        legend.position = "right") +
  labs(title = "Batch01 DAG - lipid species concentrations per population", y = "Log concentrations", x = "Lipid species (Number of carbon _ Number of double bonds)") + 
  theme(axis.text = element_text(size = 8, face = "bold"), 
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10), 
        legend.text = element_text(size = 11), 
        legend.title = element_text(size = 11))  


batch01.plot.DG
ggsave("batch01.plot.DG.png", width =10, height = 8)

#pqn normalised
batch01.plot.DG.pqn <-  ggplot(batch01.pop %>% filter(Class == "DG"), aes(paste(`Lipid species`, Name, sep = " "), log(pqn), fill = `population`)) +
  geom_boxplot(width =1) +
  scale_fill_brewer(palette = "Set1") +
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text = element_text(size = 7, vjust = -0.25), 
        axis.text.x = element_text(angle = 40), 
        panel.background = element_rect(fill = "grey91"), 
        legend.position = "right") +
  labs(title = "Batch01 DAG - lipid species normalised area per population", y = "Log(pqn)", x = "Lipid species (Number of carbon _ Number of double bonds)") + 
  theme(axis.text = element_text(size = 8, face = "bold"), 
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10), 
        legend.text = element_text(size = 11), 
        legend.title = element_text(size = 11))  


batch01.plot.DG.pqn

##%######################################################%##
#                                                          #
####                      PS plot                       ####
#                                                          #
##%######################################################%##

batch01.PS <-  batch01.pop %>% 
  filter(Class == "PS") 

batch01.plot.PS <-  ggplot(batch01.PS, aes(`Lipid species`, log(conc_mgmL_compounds), fill = `population`)) +
  geom_boxplot(width =1) +
  scale_fill_brewer(palette = "Set1") +
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text = element_text(size = 7, vjust = -0.25), 
        axis.text.x = element_text(angle = 40), 
        panel.background = element_rect(fill = "grey91"), 
        legend.position = "right") +
  labs(title = "Batch01 PS - lipid species concentrations per population", y = "Log concentrations", x = "Lipid species (Number of carbon _ Number of double bonds)") + 
  theme(axis.text = element_text(size = 8, face = "bold"), 
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10), 
        legend.text = element_text(size = 11), 
        legend.title = element_text(size = 11))  


batch01.plot.PS
ggsave("batch01.plot.PS.png", width =10, height = 8)

#pqn normalised
batch01.plot.PS.pqn <-  ggplot(batch01.PS, aes(`Lipid species`, log(pqn), fill = `population`)) +
  geom_boxplot(width =1) +
  scale_fill_brewer(palette = "Set1") +
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text = element_text(size = 7, vjust = -0.25), 
        axis.text.x = element_text(angle = 40), 
        panel.background = element_rect(fill = "grey91"), 
        legend.position = "right") +
  labs(title = "Batch01 PS - lipid species normalised area per population", y = "Log(pqn)", x = "Lipid species (Number of carbon _ Number of double bonds)") + 
  theme(axis.text = element_text(size = 8, face = "bold"), 
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10), 
        legend.text = element_text(size = 11), 
        legend.title = element_text(size = 11))  


batch01.plot.PS.pqn

##%######################################################%##
#                                                          #
####                      PI plot                       ####
#                                                          #
##%######################################################%##

batch01.PI <-  batch01.pop %>% 
  filter(Class == "PI") 

batch01.plot.PI <-  ggplot(batch01.PI, aes(`Lipid species`, log(conc_mgmL_compounds), fill = `population`)) +
  geom_boxplot(width =1) +
  scale_fill_brewer(palette = "Set1") +
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text = element_text(size = 7, vjust = -0.25), 
        axis.text.x = element_text(angle = 40), 
        panel.background = element_rect(fill = "grey91"), 
        legend.position = "right") +
  labs(title = "Batch01 PI - lipid species concentrations per population", y = "Log concentrations", x = "Lipid species (Number of carbon _ Number of double bonds)") + 
  theme(axis.text = element_text(size = 8, face = "bold"), 
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10), 
        legend.text = element_text(size = 11), 
        legend.title = element_text(size = 11))  


batch01.plot.PI
ggsave("batch01.plot.PI.png", width =10, height = 8)

#pqn normalised
batch01.plot.PI.pqn <-  ggplot(batch01.PI, aes(`Lipid species`, log(pqn), fill = `population`)) +
  geom_boxplot(width =1) +
  scale_fill_brewer(palette = "Set1") +
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text = element_text(size = 7, vjust = -0.25), 
        axis.text.x = element_text(angle = 40), 
        panel.background = element_rect(fill = "grey91"), 
        legend.position = "right") +
  labs(title = "Batch01 PI - lipid species normalised area per population", y = "Log(pqn)", x = "Lipid species (Number of carbon _ Number of double bonds)") + 
  theme(axis.text = element_text(size = 8, face = "bold"), 
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10), 
        legend.text = element_text(size = 11), 
        legend.title = element_text(size = 11))  


batch01.plot.PI.pqn

##%######################################################%##
#                                                          #
####                      PG plots                      ####
#                                                          #
##%######################################################%##

batch01.PG <-  batch01.pop %>% 
  filter(Class == "PG") 

batch01.plot.PG <-  ggplot(batch01.PG, aes(`Lipid species`, log(conc_mgmL_compounds), fill = `population`)) +
  geom_boxplot(width =1) +
  scale_fill_brewer(palette = "Set1") +
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text = element_text(size = 7, vjust = -0.25), 
        axis.text.x = element_text(angle = 40), 
        panel.background = element_rect(fill = "grey91"), 
        legend.position = "right") +
  labs(title = "Batch01 PG - lipid species concentrations per population", y = "Log concentrations", x = "Lipid species (Number of carbon _ Number of double bonds)") + 
  theme(axis.text = element_text(size = 8, face = "bold"), 
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10), 
        legend.text = element_text(size = 11), 
        legend.title = element_text(size = 11))  


batch01.plot.PG
ggsave("batch01.plot.PG.png", width =10, height = 8)

#pqn normalised
batch01.plot.PG.pqn <-  ggplot(batch01.PG, aes(`Lipid species`, log(pqn), fill = `population`)) +
  geom_boxplot(width =1) +
  scale_fill_brewer(palette = "Set1") +
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text = element_text(size = 7, vjust = -0.25), 
        axis.text.x = element_text(angle = 40), 
        panel.background = element_rect(fill = "grey91"), 
        legend.position = "right") +
  labs(title = "Batch01 PG - lipid species normalised area per population", y = "Log(pqn)", x = "Lipid species (Number of carbon _ Number of double bonds)") + 
  theme(axis.text = element_text(size = 8, face = "bold"), 
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10), 
        legend.text = element_text(size = 11), 
        legend.title = element_text(size = 11))  


batch01.plot.PG.pqn

##%######################################################%##
#                                                          #
####                      TAG plot                       ####
#                                                          #
##%######################################################%##

#sum the TAG species
batch01.TG.sum <-  batch01.pop %>% 
  filter(Class == "TG") %>% 
  group_by(Class, `Lipid species`, population, sample) %>% 
  summarise(sumpqn = sum(pqn),
         sum_conc_mgml = sum(conc_mgmL_compounds))

batch01.plot.TG <-  ggplot(batch01.TG.sum, aes(`Lipid species`, log(sum_conc_mgml), fill = `population`)) +
  geom_boxplot(width =1) +
  scale_fill_brewer(palette = "Set1") +
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text = element_text(size = 7, vjust = -0.25), 
        axis.text.x = element_text(angle = 40), 
        panel.background = element_rect(fill = "grey91"), 
        legend.position = "right") +
  labs(title = "Batch01 TAG - lipid species concentrations per population", y = "Log concentrations", x = "Lipid species (Number of carbon _ Number of double bonds)") + 
  theme(axis.text = element_text(size = 8, face = "bold"), 
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10), 
        legend.text = element_text(size = 11), 
        legend.title = element_text(size = 11)) 
 
batch01.plot.TG
ggsave("batch01.plot.TG.png", width =15, height = 8)

#pqn normalisation
batch01.plot.TG.pqn <-  ggplot(batch01.TG.sum, aes(`Lipid species`, log(sumpqn), fill = `population`)) +
  geom_boxplot(width =1) +
  scale_fill_brewer(palette = "Set1") +
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text = element_text(size = 7, vjust = -0.25), 
        axis.text.x = element_text(angle = 40), 
        panel.background = element_rect(fill = "grey91"), 
        legend.position = "right") +
  labs(title = "Batch01 TAG - lipid species normalised area per population", y = "Log(pqn)", x = "Lipid species (Number of carbon _ Number of double bonds)") + 
  theme(axis.text = element_text(size = 8, face = "bold"), 
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10), 
        legend.text = element_text(size = 11), 
        legend.title = element_text(size = 11)) 

batch01.plot.TG.pqn

##%######################################################%##
#                                                          #
####                      PC plot                       ####
#                                                          #
##%######################################################%##

batch01.PC <-  batch01.pop %>% 
  filter(Class == "PC") 

batch01.plot.PC <-  ggplot(batch01.PC, aes(`Lipid species`, log(conc_mgmL_compounds), fill = `population`)) +
  geom_boxplot(width =1) +
  scale_fill_brewer(palette = "Set1") +
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text = element_text(size = 7, vjust = -0.25), 
        axis.text.x = element_text(angle = 40), 
        panel.background = element_rect(fill = "grey91"), 
        legend.position = "right") +
  labs(title = "Batch01 PC - lipid species concentrations per population", y = "Log concentrations", x = "Lipid species (Number of carbon _ Number of double bonds)") + 
  theme(axis.text = element_text(size = 8, face = "bold"), 
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10), 
        legend.text = element_text(size = 11), 
        legend.title = element_text(size = 11))  


batch01.plot.PC
ggsave("batch01.plot.PC.png", width =10, height = 8)


#pqn normalisation
batch01.plot.PC.pqn <-  ggplot(batch01.PC, aes(`Lipid species`, log(pqn), fill = `population`)) +
  geom_boxplot(width =1) +
  scale_fill_brewer(palette = "Set1") +
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.text = element_text(size = 7, vjust = -0.25), 
        axis.text.x = element_text(angle = 40), 
        panel.background = element_rect(fill = "grey91"), 
        legend.position = "right") +
  labs(title = "Batch01 PC - lipid species normalised area per population", y = "Log(pqn)", x = "Lipid species (Number of carbon _ Number of double bonds)") + 
  theme(axis.text = element_text(size = 8, face = "bold"), 
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10), 
        legend.text = element_text(size = 11), 
        legend.title = element_text(size = 11))  


batch01.plot.PC.pqn
# KABLE ###############################
batch01.table <- batch01.values %>% 
  kable() %>% 
  kable_styling()

batch01.table
#######################################

# Substituting NA in Oleic Acid: Lipid species, TC, DB (in p1)

#ylab("Number of Flies Died")+xlab("Number of hours after initiation of desiccation stress")+ylim(0,80)+#ylim and xlim fxn changes x & y limits
#scale_x_continuous(breaks = seq(0,60,by=5))#sets spacing in x-axis

# test01 <- p1 %>% 
#   select (TC) 



 # filter(Name == "Oleic acid")

#is.na(x)) %>% 
  %>% 
 # grepl (TC =="NA", "18", " ")

#  test01

