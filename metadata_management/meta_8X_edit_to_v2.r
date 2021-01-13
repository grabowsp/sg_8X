# Addition of UNI_ACC column to lump together samples from different collectors
#   but should have originated from the same source population

# on HA
# bash
# source activate R_analysis

### LOAD PACKAGES ###
library(data.table)

### LOAD DATA ###
meta_file <- '/home/f2p1/work/grabowsk/data/switchgrass/sg_8X_metadata_v1.0.csv'
samp_meta <- fread(meta_file)

### SET OUTPUT ###

out_file <- '/home/f2p1/work/grabowsk/data/switchgrass/sg_8X_metadata_v2.0.csv'

#################
weird_locality_inds <- c(53,797,921,924,936)

samp_meta[ , UNI_ACC := as.character(NA)]


look_for_acc <- function(search_string){
  si1 <- grep(search_string, samp_meta$SOURCE_ID_1)
  si2 <- grep(search_string, samp_meta$SOURCE_ID_2)
  l1 <- grep(search_string, samp_meta$LOCALITY)
  return(unique(c(si1, si2, l1)))
}

samp_meta[ACC == sort(unique(samp_meta$ACC[which(is.na(samp_meta$UNI_ACC))]))[1], ]

samp_meta[ACC == 'J001', UNI_ACC := 'PI_422016']
look_for_acc('422016')

samp_meta[ACC == 'J002', UNI_ACC := 'CITRUS_CO_FL']
look_for_acc('Citrus')

samp_meta[ACC == 'J003', ]
samp_meta[look_for_acc('315723'), UNI_ACC := 'HOFFMAN']

samp_meta[ACC == 'J004', UNI_ACC := 'ELLSWORTH']
look_for_acc('315724')
look_for_acc('BN-10860-61')
look_for_acc('Ellsworth')

samp_meta[ACC == 'J005', ]
samp_meta[look_for_acc('14669'), UNI_ACC := 'PI_315725']
look_for_acc('315725')
look_for_acc('BN-14669-92')
look_for_acc('14669')
look_for_acc('Coffeeville')

samp_meta[ACC == 'J006', UNI_ACC := 'PI_315727']
look_for_acc('315727')
look_for_acc('11357')

samp_meta[ACC == 'J007', UNI_ACC := 'ARGENTINA']
look_for_acc('337553')
look_for_acc('Argentina')
look_for_acc('Rafaela')

samp_meta[ACC == 'J008', UNI_ACC := 'PI_414065']
look_for_acc('414065')
look_for_acc('14668')
#look_for_acc('Pangburn')

samp_meta[ACC == 'J009', UNI_ACC := 'PI_414066']
look_for_acc('414066')
look_for_acc('5669')
look_for_acc('Grenville')

samp_meta[ACC == 'J010', UNI_ACC := 'PI_414067']
look_for_acc('414067')

samp_meta[ACC == sort(unique(samp_meta$ACC[which(is.na(samp_meta$UNI_ACC))]))[1], ]
#samp_meta[ACC == 'J012', UNI_ACC := ]

samp_meta[ACC == 'J011', UNI_ACC := 'PI_414068']
look_for_acc('414068')
look_for_acc('18758')

samp_meta[ACC == 'J012', UNI_ACC := 'PI_414069']
look_for_acc('414069')
look_for_acc('309')

samp_meta[ACC == 'J013', UNI_ACC := 'PI_414070']
look_for_acc('414070')
look_for_acc('12323')

samp_meta[ACC == 'J014', UNI_ACC := 'CARTHAGE']
samp_meta[ACC == 'J438', UNI_ACC := 'CARTHAGE']
look_for_acc('421138')
look_for_acc('Carthage')
look_for_acc('CARTHAGE')

samp_meta[ACC == 'J015', UNI_ACC := 'PI_421520']
look_for_acc('421520')

samp_meta[ACC == 'J016', UNI_ACC := ]
samp_meta[look_for_acc('Kanlow'), UNI_ACC := 'KANLOW']
look_for_acc('Kanlow')
samp_meta[VCF_NAME == 'NL94', UNI_ACC := as.character(NA)]
# VCF_NAME == 'NL94' is NOT exactly Kanlow
#look_for_acc('KANLOW')
#look_for_acc('421521')
#look_for_acc('Wetumka')
#look_for_acc('KAN')

samp_meta[ACC == 'J017', UNI_ACC := 'MIAMI']
look_for_acc('421901')
look_for_acc('Miami')
look_for_acc('MIAMI')

samp_meta[ACC == 'J018', UNI_ACC := 'PI_421999']
look_for_acc('42199')
look_for_acc('AM-314')
#look_for_acc('314')

samp_meta[ACC == 'J019', UNI_ACC := 'PI_422000']
look_for_acc('422000')
look_for_acc('F-686')
look_for_acc('686')
look_for_acc('Wabasso')

samp_meta[ACC == 'J020', UNI_ACC := ]
# J020 is the weird FL upland 4X; I think it's maybe supposed to be "Stuart"
#  but certainly isn't - need to check this
samp_meta[look_for_acc('422001'), UNI_ACC := 'PI_422001']
#look_for_acc('Stuart')
#look_for_acc('STUART')

samp_meta[ACC == sort(unique(samp_meta$ACC[which(is.na(samp_meta$UNI_ACC))]))[1], ]
#samp_meta[ACC == 'J012', UNI_ACC := ]

samp_meta[ACC == 'J021', UNI_ACC := 'PI_422003']
look_for_acc('422003')
look_for_acc('PMT-785')
look_for_acc('785')

samp_meta[ACC == 'J022', UNI_ACC := ]
look_for_acc('422006')
look_for_acc('Alamo')
look_for_acc('ALAMO')
samp_meta[unique(c(look_for_acc('422006'), look_for_acc('Alamo'), 
  look_for_acc('ALAMO'))), UNI_ACC := 'ALAMO']
samp_meta[ACC == 'J287', UNI_ACC := 'ALAMO']
# this includes a few AP13 samples: J287_A, J206.A, J223.A

samp_meta[ACC == 'J023', UNI_ACC := 'KY_1625']
look_for_acc('431575')
look_for_acc('1625')

samp_meta[ACC == 'J025', UNI_ACC := 'CAVE_IN_ROCK']
samp_meta[ACC == 'J439', UNI_ACC := 'CAVE_IN_ROCK']
samp_meta[VCF_NAME == 'Cave_In_Rock_WO1', UNI_ACC := 'CAVE_IN_ROCK']
look_for_acc('469228')
look_for_acc('Cave')
look_for_acc('Rock')
# J025, J439, VCF_NAME == Cave_In_Rock_WO1

samp_meta[ACC == 'J026', UNI_ACC := 'PI_476290']
look_for_acc('476290')
look_for_acc('2086')

samp_meta[ACC == 'J027', UNI_ACC := 'PI_476291']
samp_meta[ACC == 'J027', COLLECTION_TYPE := 'Breeding Selection']
look_for_acc('476291')
look_for_acc('2099')

samp_meta[ACC == 'J028', UNI_ACC := 'PI_476292']
look_for_acc('476292')
look_for_acc(2100)

samp_meta[ACC == 'J029', UNI_ACC := 'PI_476293']
look_for_acc('476293')
look_for_acc('2101')

load(tmp_workspace)

samp_meta[ACC == 'J030', UNI_ACC := ]
# turns out this is two separate GRIN accessions... continue from here
samp_meta[look_for_acc('476295'), UNI_ACC := 'PI_476295']
samp_meta[look_for_acc('476295'), LOCALITY := 'Colorado Springs, Colorado']
samp_meta[UNI_ACC == 'PI_476295', LATITUDE := 38.92]
samp_meta[UNI_ACC == 'PI_476295', LONGITUDE := -104.67]
samp_meta[UNI_ACC == 'PI_476295', NOTE_LATLONG := 'Generic Colorado Springs coordinates']

samp_meta[look_for_acc('476294'), UNI_ACC := 'PI_476294']

samp_meta[ACC == sort(unique(samp_meta$ACC[which(is.na(samp_meta$UNI_ACC))]))[1], ]
#samp_meta[ACC == 'J012', UNI_ACC := ]

samp_meta[ACC == 'J032', UNI_ACC := 'PI_476296']
samp_meta[ACC == 'J032', LATITUDE := 38.78]
samp_meta[ACC == 'J032', LONGITUDE := -104.71]
samp_meta[ACC == 'J032', NOTE_LATLONG := 'generic Colorado Springs coordinates']
samp_meta[ACC == 'J032', STATE := 'Colorado']
samp_meta[ACC == 'J032', LOCALITY := 'Colorado Springs, Colorado']
samp_meta[ACC == 'J032', ELEVATION := 1839]

look_for_acc('476296')

samp_meta[ACC == 'J033', UNI_ACC := 'PI_476297']
look_for_acc('476297')
look_for_acc('Caddo')

samp_meta[ACC == 'J036', UNI_ACC := 'PI_478002']
look_for_acc('478002')
look_for_acc('6011')

samp_meta[ACC == 'J037', UNI_ACC := ]
look_for_acc('537588')
samp_meta[look_for_acc('Dacotah'), UNI_ACC := 'DACOTAH']
samp_meta[VCF_NAME == 'DAC6' , UNI_ACC := 'DACOTAH']
samp_meta[UNI_ACC == 'DACOTAH', COLLECTION_TYPE := 'Cultivar']

samp_meta[ACC == 'J038', UNI_ACC := 'TRAILBLAZER']
look_for_acc('549094')
look_for_acc('TrailBlazer')

samp_meta[ACC == 'J039', UNI_ACC := 'SHAWNEE']
look_for_acc('591824')
look_for_acc('SHAWNEE')
look_for_acc('Shawnee')

samp_meta[ACC == 'J040', ]
look_for_acc('598136')
samp_meta[look_for_acc('Sunburst'), UNI_ACC := 'SUNBURST']

samp_meta[ACC == 'J041', UNI_ACC := 'TEM-LODORM']
look_for_acc('636468')
look_for_acc('LoDorm')
look_for_acc('TEM')

samp_meta[ACC == 'J043', UNI_ACC := 'WS8U']
look_for_acc('639192')
look_for_acc('WS8U')

samp_meta[ACC == 'J044', UNI_ACC := 'FALCON']
look_for_acc('642190')

samp_meta[ACC == 'J045', UNI_ACC := 'SUMMER']
look_for_acc('642191')
samp_meta[look_for_acc('Summer'), UNI_ACC := 'SUMMER']
samp_meta[UNI_ACC == 'SUMMER', COLLECTION_TYPE := 'Cultivar']

samp_meta[ACC == sort(unique(samp_meta$ACC[which(is.na(samp_meta$UNI_ACC))]))[1], ]
#samp_meta[ACC == 'J012', UNI_ACC := ]

samp_meta[ACC == 'J046', ]
look_for_acc('642192')
samp_meta[look_for_acc('Pathfinder'), UNI_ACC := 'PATHFINDER']

samp_meta[ACC == 'J058', UNI_ACC := 'ND_ARVID_BOE']
look_for_acc('642207')
samp_meta[ACC == 'J066', UNI_ACC := 'ND_ARVID_BOE']
look_for_acc('642218')
samp_meta[ACC == 'J067', UNI_ACC := 'ND_ARVID_BOE']
look_for_acc('642219')
samp_meta[ACC == 'J068', UNI_ACC := 'ND_ARVID_BOE']

samp_meta[ACC == 'J069', UNI_ACC := 'OSSP_FL']
look_for_acc('OSSP')
look_for_acc('Scherer')

samp_meta[ACC == 'J070', UNI_ACC := 'ND_ARVID_BOE']
samp_meta[ACC == 'J071', UNI_ACC := 'ND_ARVID_BOE']

samp_meta[ACC == 'J072', UNI_ACC := 'PASCO_FL']
look_for_acc('Pasco')

samp_meta[ACC == 'J076', UNI_ACC := 'ND_ARVID_BOE']
samp_meta[ACC == 'J077', UNI_ACC := 'ND_ARVID_BOE']
samp_meta[ACC == 'J078', UNI_ACC := 'ND_ARVID_BOE']
samp_meta[ACC == 'J079', UNI_ACC := 'ND_ARVID_BOE']
samp_meta[ACC == 'J080', UNI_ACC := 'ND_ARVID_BOE']
samp_meta[ACC == 'J081', UNI_ACC := 'ND_ARVID_BOE']
samp_meta[ACC == 'J082', UNI_ACC := 'ND_ARVID_BOE']

samp_meta[ACC == 'J083', UNI_ACC := 'SP_BLUFF_1_2']
look_for_acc('SP Bluff')
look_for_acc('Sprewell')

samp_meta[ACC == 'J097', UNI_ACC := 'ND_ARVID_BOE']
samp_meta[ACC == 'J098', UNI_ACC := 'ND_ARVID_BOE']
samp_meta[ACC == 'J103', UNI_ACC := 'ND_ARVID_BOE']
samp_meta[ACC == 'J104', UNI_ACC := 'ND_ARVID_BOE']
samp_meta[ACC == 'J105', UNI_ACC := 'ND_ARVID_BOE']
samp_meta[ACC == 'J107', UNI_ACC := 'ND_ARVID_BOE']
samp_meta[ACC == 'J109', UNI_ACC := 'ND_ARVID_BOE']
samp_meta[ACC == 'J110', UNI_ACC := 'ND_ARVID_BOE_71SG']

samp_meta[UNI_ACC == 'ND_ARVID_BOE', UNI_ACC := 'ND_ARVID_BOE_70SG']

samp_meta[ACC == sort(unique(samp_meta$ACC[which(is.na(samp_meta$UNI_ACC))]))[1], ]
#samp_meta[ACC == 'J012', UNI_ACC := ]

samp_meta[ACC == 'J111', UNI_ACC := 'ND_ARVID_BOE_71SG']
samp_meta[ACC == 'J112', UNI_ACC := 'ND_ARVID_BOE_71SG']
samp_meta[ACC == 'J113', UNI_ACC := 'ND_ARVID_BOE_71SG']
samp_meta[ACC == 'J115', UNI_ACC := 'ND_ARVID_BOE_71SG']
samp_meta[ACC == 'J116', UNI_ACC := 'ND_ARVID_BOE_71SG']
samp_meta[ACC == 'J117', UNI_ACC := 'ND_ARVID_BOE_71SG']
samp_meta[ACC == 'J118', UNI_ACC := 'ND_ARVID_BOE_71SG']
samp_meta[ACC == 'J119', UNI_ACC := 'ND_ARVID_BOE_71SG']
samp_meta[grep('71SG', samp_meta$SOURCE_ID_2), UNI_ACC := 'ND_ARVID_BOE_71SG']
samp_meta[grep('70SG', samp_meta$SOURCE_ID_2), UNI_ACC := 'ND_ARVID_BOE_70SG']

samp_meta[look_for_acc('Blackwell'), UNI_ACC := 'BLACKWELL']
look_for_acc('657661')
look_for_acc('657663')
look_for_acc('Blackwell')

samp_meta[ACC == 'J164', UNI_ACC := 'PI_659333']
look_for_acc('659333')

samp_meta[ACC == 'J165', UNI_ACC := 'PI_659334']
look_for_acc('659334')
look_for_acc('2008FL')

samp_meta[ACC == 'J166', UNI_ACC := 'PI_659335']
look_for_acc('659335')

samp_meta[ACC == 'J169', UNI_ACC := 'PI_659340']
look_for_acc('659340')
look_for_acc('9064231')

samp_meta[ACC == 'J170', UNI_ACC := 'PI_659341']
look_for_acc('659341')

samp_meta[ACC == 'J171', UNI_ACC := 'PI_659342']
look_for_acc('659342')

samp_meta[ACC == 'J172', UNI_ACC := 'PI_659343']
look_for_acc('659343')

samp_meta[ACC == 'J173', UNI_ACC := 'PI_659344']
look_for_acc('659344')

samp_meta[ACC == 'J174', UNI_ACC := 'PI_659345']
look_for_acc('659345')

samp_meta[ACC == sort(unique(samp_meta$ACC[which(is.na(samp_meta$UNI_ACC))]))[1], ]
#samp_meta[ACC == 'J012', UNI_ACC := ]

samp_meta[ACC == 'J175', UNI_ACC := 'PI_659346']
look_for_acc('659346')

samp_meta[ACC == 'J177', UNI_ACC := 'PI_666207']
look_for_acc('666207')

samp_meta[ACC == 'J178', UNI_ACC := 'PI_666208']
look_for_acc('666208')

samp_meta[ACC == 'J179', UNI_ACC := 'PI_666209']
look_for_acc('666209')

samp_meta[ACC == 'J180', UNI_ACC := 'PI_666210']
look_for_acc('666210')

samp_meta[ACC == 'J181', UNI_ACC := 'PI_666211']
look_for_acc('666211')

samp_meta[ACC == 'J182', UNI_ACC := 'PI_666212']
look_for_acc('666212')

samp_meta[ACC == 'J184', UNI_ACC := 'PI_666214']
look_for_acc('666214')

samp_meta[ACC == 'J185', UNI_ACC := 'PI_666215']
look_for_acc('666215')

samp_meta[ACC == 'J186', UNI_ACC := 'PI_666216']
look_for_acc('666216')

samp_meta[ACC == 'J187', UNI_ACC := 'PI_666217']
look_for_acc('666217')

samp_meta[ACC == 'J188', UNI_ACC := 'PI_666218']
look_for_acc('666218')

samp_meta[ACC == 'J189', UNI_ACC := 'PI_666219']
look_for_acc('666219')

samp_meta[ACC == 'J190', UNI_ACC := 'PI_666220']
look_for_acc('666220')

samp_meta[ACC == 'J191', UNI_ACC := 'PI_666221']
look_for_acc('666221')

samp_meta[ACC == 'J192', UNI_ACC := 'PI_667580']
look_for_acc('667580')

samp_meta[ACC == sort(unique(samp_meta$ACC[which(is.na(samp_meta$UNI_ACC))]))[1], ]
#samp_meta[ACC == 'J012', UNI_ACC := ]

samp_meta[ACC == 'J193', UNI_ACC := 'PI_667581']
look_for_acc('667581')

samp_meta[ACC == 'J194', UNI_ACC := 'PI_667582']
look_for_acc('667582')

samp_meta[ACC == 'J195', UNI_ACC := 'PI_667583']
look_for_acc('667583')

samp_meta[ACC == 'J197', UNI_ACC := 'PI_667585']
look_for_acc('667585')

samp_meta[ACC == 'J198', UNI_ACC := 'PI_667586']
look_for_acc('667586')

samp_meta[ACC == 'J199', UNI_ACC := 'PI_667587']
look_for_acc('667587')

samp_meta[ACC == 'J200', UNI_ACC := 'PI_667588']
look_for_acc('667588')

samp_meta[ACC == 'J201', UNI_ACC := 'PI_667589']
look_for_acc('667589')

samp_meta[ACC == 'J202', UNI_ACC := 'PI_667590']
look_for_acc('667590')

samp_meta[ACC == 'J204', UNI_ACC := 'PI_667592']
look_for_acc('667592')

samp_meta[ACC == 'J205', UNI_ACC := 'PI_667593']
look_for_acc('667593')

samp_meta[ACC == 'J208', UNI_ACC := 'PI_315728']
look_for_acc('315728')
look_for_acc('13645')

samp_meta[ACC == 'J209', UNI_ACC := 'FUL_1']
look_for_acc('FUL')
samp_meta[ACC == 'J294', UNI_ACC := 'FUL_2']

samp_meta[ACC == 'J210', UNI_ACC := 'DVR']
look_for_acc('DVR')

samp_meta[ACC == sort(unique(samp_meta$ACC[which(is.na(samp_meta$UNI_ACC))]))[1], ]
#samp_meta[ACC == 'J012', UNI_ACC := ]

samp_meta[ACC == 'J211', UNI_ACC := 'TOB']
look_for_acc('TOB')

samp_meta[ACC == 'J212', UNI_ACC := 'MIS']
look_for_acc('MIS')

samp_meta[ACC == 'J213', UNI_ACC := 'ENC']
look_for_acc('ENC')

samp_meta[ACC == 'J214', UNI_ACC := 'NAS']
look_for_acc('NAS')

#samp_meta[ACC == 'J215', UNI_ACC := ]
samp_meta[look_for_acc('9090293'), UNI_ACC := 'FAN']
look_for_acc('FAN')

samp_meta[ACC == 'J216', UNI_ACC := 'HCC']
look_for_acc('HCC')

samp_meta[ACC == 'J217', UNI_ACC := 'WAR']
look_for_acc('WAR')

samp_meta[ACC == 'J218', UNI_ACC := 'OAS']
look_for_acc('OAS')

samp_meta[ACC == 'J219', UNI_ACC := 'WIL']
look_for_acc('WIL')

samp_meta[ACC == 'J220', UNI_ACC := 'PLV']
look_for_acc('PLV')

#samp_meta[ACC == 'J221', UNI_ACC := ]
samp_meta[look_for_acc('9086192'), UNI_ACC := 'FLW']
look_for_acc('FLW')
look_for_acc('9086192')

samp_meta[ACC == 'J222', UNI_ACC := 'DEV']
look_for_acc('DEV')

samp_meta[ACC == 'J224', UNI_ACC := 'LAS']
look_for_acc('LAS')

samp_meta[ACC == 'J225', UNI_ACC := 'BRO']
look_for_acc('BRO')

samp_meta[ACC == 'J226', UNI_ACC := 'HOD']
look_for_acc('HOD')

#samp_meta[ACC == 'J228', UNI_ACC := ]
samp_meta[look_for_acc(9088695), UNI_ACC := 'ANW']
look_for_acc('ANW')

samp_meta[ACC == sort(unique(samp_meta$ACC[which(is.na(samp_meta$UNI_ACC))]))[1], ]
#samp_meta[ACC == 'J012', UNI_ACC := ]

samp_meta[ACC == 'J229', UNI_ACC := 'LMV']
look_for_acc('LMV')

#samp_meta[ACC == 'J230', UNI_ACC := ]
samp_meta[look_for_acc('WBC'), UNI_ACC := 'WBC']
look_for_acc('WBC')

samp_meta[ACC == 'J231', UNI_ACC := 'DHR']
look_for_acc('DHR')

samp_meta[ACC == 'J232', UNI_ACC := 'GRN']
look_for_acc('GRN')

samp_meta[ACC == 'J234', UNI_ACC := 'DOE']
look_for_acc('DOE')

samp_meta[ACC == 'J235', UNI_ACC := 'ATM']
look_for_acc('ATM')

#samp_meta[ACC == 'J236', UNI_ACC := ]
samp_meta[look_for_acc('9086194'), UNI_ACC := 'FLO']
look_for_acc('FLO')

samp_meta[ACC == 'J237', UNI_ACC := 'HDR']
look_for_acc('HDR')

samp_meta[ACC == 'J238', UNI_ACC := 'CMH']
look_for_acc('CMH')

samp_meta[ACC == 'J240', UNI_ACC := 'COP']
look_for_acc('COP')

samp_meta[ACC == 'J241', UNI_ACC := 'SEL']
look_for_acc('SEL')

samp_meta[ACC == 'J245', UNI_ACC := 'BRA']
look_for_acc('BRA')

samp_meta[ACC == 'J246', UNI_ACC := 'PPD']
look_for_acc('PPD')

samp_meta[ACC == 'J247', UNI_ACC := 'LIH']
look_for_acc('LIH')

samp_meta[ACC == sort(unique(samp_meta$ACC[which(is.na(samp_meta$UNI_ACC))]))[1], ]
#samp_meta[ACC == 'J012', UNI_ACC := ]

#samp_meta[ACC == 'J248', UNI_ACC := ]
samp_meta[look_for_acc('9086193'), UNI_ACC := 'ARM']
look_for_acc('ARM')

samp_meta[ACC == 'J249', UNI_ACC := 'HMP']
look_for_acc('HMP')

samp_meta[ACC == 'J250', UNI_ACC := 'RMR']
look_for_acc('RMR')

samp_meta[ACC == 'J251', UNI_ACC := 'SHC']
look_for_acc('SHC')

samp_meta[ACC == 'J252', UNI_ACC := 'SW384']
look_for_acc('SW384')
look_for_acc('Ipswitch')

samp_meta[ACC == 'J253', UNI_ACC := 'SW1983']
look_for_acc('SW1983')

samp_meta[ACC == 'J254', UNI_ACC := 'SW2121']
look_for_acc('SW2121')

samp_meta[ACC == 'J256', UNI_ACC := 'SW2122']
look_for_acc('SW2122')

samp_meta[ACC == 'J257', UNI_ACC := 'BLACK_DOG']
look_for_acc('Black Dog')

samp_meta[ACC == 'J258', UNI_ACC := 'HELEN_ALLISON']
look_for_acc('Helen')

samp_meta[ACC == 'J259', UNI_ACC := 'KASOTA']
look_for_acc('Kasota')

samp_meta[ACC == 'J262', UNI_ACC := 'LUNDBLAD']
look_for_acc('Lundblad')

samp_meta[ACC == 'J263', UNI_ACC := 'PRAIRIE_COTEAU']
look_for_acc('Coteau')

samp_meta[ACC == 'J264', UNI_ACC := 'RUSHFORD']
look_for_acc('Rushford')

samp_meta[ACC == 'J265', UNI_ACC := 'ST_CROIX']
look_for_acc('Croix')

samp_meta[ACC == 'J267', UNI_ACC := 'UNCAS']
look_for_acc('Uncas')

#samp_meta[ACC == 'J267', UNI_ACC := ]
samp_meta[look_for_acc('MEX.239'), UNI_ACC := 'MEX_239']
look_for_acc('MEX')

samp_meta[ACC == sort(unique(samp_meta$ACC[which(is.na(samp_meta$UNI_ACC))]))[1], ]
#samp_meta[ACC == 'J012', UNI_ACC := ]

#samp_meta[ACC == 'J272', UNI_ACC := ]
samp_meta[look_for_acc('MEX.244'), UNI_ACC := 'MEX_244']
look_for_acc('MEX.244')

#samp_meta[ACC == 'J275', UNI_ACC := ]
samp_meta[look_for_acc('MEX.245'), UNI_ACC := 'MEX_245']
look_for_acc('MEX.245')
samp_meta[ACC == 'J276', UNI_ACC := 'MEX_245']

samp_meta[ACC == 'J279', UNI_ACC := 'MEX_247']
look_for_acc('MEX.247')

#samp_meta[ACC == 'J280', UNI_ACC := ]
samp_meta[look_for_acc('Devil'), UNI_ACC := 'DVR']
look_for_acc('Devil')

samp_meta[ACC == 'J282', UNI_ACC := 'RIVERFRONT_PARK']
look_for_acc('Riverfront')

samp_meta[ACC == 'J283', UNI_ACC := 'MARKHAM']
look_for_acc('Markham')

samp_meta[ACC == 'J285', UNI_ACC := 'HAMMONASSET']
look_for_acc('Hammonasset')

samp_meta[ACC == 'J290', UNI_ACC := 'ALAMO']
look_for_acc('A4')

#samp_meta[ACC == 'J293', UNI_ACC := ]
samp_meta[look_for_acc('9093168'), UNI_ACC := 'BEX5']
look_for_acc('BEX')

samp_meta[ACC == 'J296', UNI_ACC := 'NAM_15_01']
look_for_acc('NAM')
look_for_acc('15_01')

samp_meta[ACC == 'J297', UNI_ACC := 'NAM_15_04']
look_for_acc('15_04')

samp_meta[ACC == 'J299', UNI_ACC := 'NAM_19_05']
look_for_acc('19_05')

samp_meta[ACC == 'J300', UNI_ACC := 'NAM_31_03']
look_for_acc('31_03')

samp_meta[ACC == 'J301', UNI_ACC := 'NAM_31_16']
look_for_acc('31_16')

samp_meta[ACC == sort(unique(samp_meta$ACC[which(is.na(samp_meta$UNI_ACC))]))[1], ]
#samp_meta[ACC == 'J012', UNI_ACC := ]

samp_meta[ACC == 'J306', UNI_ACC := 'NEAR_DHR']
# is very close to Double Helix Ranch (DHR)
look_for_acc('JP')

samp_meta[ACC == 'J307', UNI_ACC := 'WIL_2']
# no collection info, in Wilson county, but there is already WIL from same
#  county, so calling WIL_2
look_for_acc('9089249')

samp_meta[ACC == 'J309', UNI_ACC := 'SAN_PATRICIO']
look_for_acc('9090297')

samp_meta[ACC == 'J312', UNI_ACC := 'KENEDY_1']
# only info I have is from Kenedy county
look_for_acc('9093263')

samp_meta[ACC == 'J313', UNI_ACC := 'HARRIS_1']
samp_meta[ACC == 'J314', UNI_ACC := 'HARRIS_2']
samp_meta[ACC == 'J340', UNI_ACC := 'HARRIS_3']
# these are all from Harris county but different collections
look_for_acc('9093270')

#samp_meta[ACC == 'J315', UNI_ACC := ]
samp_meta[look_for_acc('GA991'), UNI_ACC := 'NAM_PARENTS_GA991']
look_for_acc('GA991')
look_for_acc('NFGA36')
samp_meta[UNI_ACC == 'NAM_PARENTS_GA991', 
  COLLECTION_TYPE := 'Breeding Selection']

#samp_meta[ACC == 'J316', UNI_ACC := ]
samp_meta[look_for_acc('GA993'), UNI_ACC := 'NAM_PARENTS_GA993']
samp_meta[UNI_ACC == 'NAM_PARENTS_GA993', 
  COLLECTION_TYPE := 'Breeding Selection']

#samp_meta[ACC == 'J317', UNI_ACC := ]
samp_meta[look_for_acc('GA992'), UNI_ACC := 'NAM_PARENTS_GA992']
look_for_acc('GA992')
samp_meta[look_for_acc('GA992'), COLLECTION_TYPE := 'Breeding Selection']

samp_meta[ACC == 'J318', UNI_ACC := 'CALHOUN_TX_1']
# in Calhoun county, TX; there are other accessions from same county
look_for_acc('9093539')

samp_meta[ACC == 'J319', UNI_ACC := 'MATAGORDA_1']
samp_meta[ACC == 'J323', UNI_ACC := 'MATAGORDA_2']
samp_meta[ACC == 'J324', UNI_ACC := 'MATAGORDA_3']
# there are several accessions from Matagorda county
look_for_acc('9093553')

samp_meta[ACC == sort(unique(samp_meta$ACC[which(is.na(samp_meta$UNI_ACC))]))[1], ]
#samp_meta[ACC == 'J012', UNI_ACC := ]

samp_meta[ACC == 'J326', UNI_ACC := 'BRAZORIA']
look_for_acc('9093628')

samp_meta[ACC == 'J327', UNI_ACC := 'VICTORIA']
look_for_acc('9109757')

samp_meta[ACC == 'J328', UNI_ACC := 'CALHOUN_TX_2']
look_for_acc('9109758')

samp_meta[ACC == 'J329', UNI_ACC := 'JACKSON_TX']
look_for_acc('9109759')

samp_meta[ACC == 'J330', UNI_ACC := 'LASALLE']
look_for_acc('9109780')

samp_meta[ACC == 'J331', UNI_ACC := 'COLORADO_TX_1']
samp_meta[ACC == 'J335', UNI_ACC := 'COLORADO_TX_2']
look_for_acc('9109790')

samp_meta[ACC == 'J336', UNI_ACC := 'AUSTIN_1']
samp_meta[ACC == 'J337', UNI_ACC := 'AUSTIN_2']
samp_meta[ACC == 'J339', UNI_ACC := 'AUSTIN_3']
look_for_acc('9109816')

samp_meta[ACC == 'J341', UNI_ACC := 'RAV']
look_for_acc('RAV')

samp_meta[ACC == 'J342', UNI_ACC := 'MAY_PRAIRIE']
look_for_acc('42493')

samp_meta[ACC == 'J343', UNI_ACC := 'KING_RANCH']
look_for_acc('9111965')

samp_meta[ACC == 'J344', UNI_ACC := 'HSPNP']
look_for_acc('HSPNP')

#samp_meta[ACC == 'J346', UNI_ACC := ]
samp_meta[look_for_acc('PBRS'), UNI_ACC := 'PBRS']
look_for_acc('PBRS')

#samp_meta[ACC == 'J348', UNI_ACC := ]
samp_meta[look_for_acc('CRP'), UNI_ACC := 'CRP']
look_for_acc('CRP')

#samp_meta[ACC == 'J350', UNI_ACC := ]
samp_meta[look_for_acc('DBOS'), UNI_ACC := 'DBOS']
look_for_acc('DBOS')

samp_meta[ACC == sort(unique(samp_meta$ACC[which(is.na(samp_meta$UNI_ACC))]))[1], ]
#samp_meta[ACC == 'J012', UNI_ACC := ]

samp_meta[ACC == 'J352', UNI_ACC := 'NERP']
look_for_acc('NERP')

samp_meta[ACC == 'J353', UNI_ACC := 'BRRP']
look_for_acc('BRRP')

samp_meta[ACC == 'J354', UNI_ACC := 'GLPSNA']
look_for_acc('GLPSNA')

#samp_meta[ACC == 'J355', UNI_ACC := ]
samp_meta[look_for_acc('MNTP'), UNI_ACC := 'MNTP']
look_for_acc('MNTP')

samp_meta[ACC == 'J358', UNI_ACC := 'DPSCA']
look_for_acc('DPSCA')

samp_meta[ACC == 'J359', UNI_ACC := 'STP']
look_for_acc('STP')

samp_meta[ACC == 'J363', UNI_ACC := 'ILP']
look_for_acc('ILP')

samp_meta[ACC == 'J365', UNI_ACC := 'TKP']
look_for_acc('TKP')

samp_meta[ACC == 'J368', UNI_ACC := 'ASGA']
look_for_acc('ASGA')

samp_meta[ACC == 'J369', UNI_ACC := 'SHIAWASSEE']
look_for_acc('RL')

samp_meta[ACC == 'J371', UNI_ACC := 'MASON_MI_1']
samp_meta[ACC == 'J380', UNI_ACC := 'MASON_MI_2']
# Only info is that it comes from Mason Co, MI; there's another accession
#  from same county

samp_meta[ACC == 'J374', UNI_ACC := 'ALLEGAN_1']
look_for_acc('AS1')

samp_meta[ACC == 'J376', UNI_ACC := 'GRAND_HAVEN']
look_for_acc('GH')

samp_meta[ACC == 'J378', UNI_ACC := 'OTTAWA_MI_1']
# not much info about the accession
look_for_acc('RM')

samp_meta[ACC == sort(unique(samp_meta$ACC[which(is.na(samp_meta$UNI_ACC))]))[1], ]
#samp_meta[ACC == 'J012', UNI_ACC := ]

samp_meta[ACC == 'J379', UNI_ACC := 'ALGONAC']
look_for_acc('AL')

samp_meta[ACC == 'J381', UNI_ACC := 'VAN_BUREN_SP']
look_for_acc('Buren')
samp_meta[ACC == 'J382', UNI_ACC := 'VAN_BUREN_AMTRAK']

samp_meta[ACC == 'J383', UNI_ACC := 'PEL']
look_for_acc('PEL')

samp_meta[ACC == 'J384', UNI_ACC := 'HLP']
look_for_acc('HLP')

samp_meta[ACC == 'J385', UNI_ACC := 'DRP']
look_for_acc('DRP')

samp_meta[ACC == 'J386', UNI_ACC := 'PBP']
look_for_acc('PBP')

samp_meta[ACC == 'J387', UNI_ACC := 'FRP']
look_for_acc('FRP')

samp_meta[ACC == 'J389', UNI_ACC := 'BHC']
look_for_acc('BHC')

samp_meta[ACC == 'J390', UNI_ACC := 'SW31']
look_for_acc('SW31')

samp_meta[ACC == 'J393', UNI_ACC := 'SW40']
look_for_acc('SW40')

samp_meta[ACC == 'J394', UNI_ACC := 'SW43']
look_for_acc('SW43')

samp_meta[ACC == 'J395', UNI_ACC := 'SW46']
look_for_acc('SW46')

samp_meta[ACC == 'J396', UNI_ACC := 'SW49']
look_for_acc('SW49')

samp_meta[ACC == 'J397', UNI_ACC := 'SW50']
look_for_acc('SW50')

samp_meta[ACC == 'J398', UNI_ACC := 'SW51']
look_for_acc('SW51')

samp_meta[ACC == sort(unique(samp_meta$ACC[which(is.na(samp_meta$UNI_ACC))]))[1], ]
#samp_meta[ACC == 'J012', UNI_ACC := ]

samp_meta[ACC == 'J399', UNI_ACC := 'SW58']
look_for_acc('SW58')

samp_meta[ACC == 'J400', UNI_ACC := 'SW63']
look_for_acc('SW63')

samp_meta[ACC == 'J402', UNI_ACC := 'SW65']
look_for_acc('SW65')

samp_meta[ACC == 'J403', UNI_ACC := 'SW102']
look_for_acc('SW102')

samp_meta[ACC == 'J406', UNI_ACC := 'SW112']
look_for_acc('SW112')

samp_meta[ACC == 'J407', UNI_ACC := 'SW114']
look_for_acc('SW114')

samp_meta[ACC == 'J409', UNI_ACC := 'SW116']
look_for_acc('SW116')

samp_meta[ACC == 'J410', UNI_ACC := 'SW122']
look_for_acc('SW122')

samp_meta[ACC == 'J411', UNI_ACC := 'SW123']
look_for_acc('SW123')

samp_meta[ACC == 'J412', UNI_ACC := 'SW124']
look_for_acc('SW124')

samp_meta[ACC == 'J413', UNI_ACC := 'SW127']
look_for_acc('SW127')

samp_meta[ACC == 'J415', UNI_ACC := 'SW129']
look_for_acc('SW129')

samp_meta[ACC == 'J416', UNI_ACC := 'SW781']
# This seems like it might be PI_674672, but not sure
look_for_acc('SW781')

samp_meta[ACC == 'J417', UNI_ACC := 'SW782']
look_for_acc('SW782')

samp_meta[ACC == 'J418', UNI_ACC := 'HIGH_TIDE']
look_for_acc('High')

samp_meta[ACC == 'J419', UNI_ACC := 'TIMBER']
samp_meta[ACC == 'J419', COLLECTION_TYPE := 'Cultivar']
look_for_acc('Timber')

samp_meta[ACC == 'J420', UNI_ACC := 'SW786']
look_for_acc('SW786')

samp_meta[ACC == sort(unique(samp_meta$ACC[which(is.na(samp_meta$UNI_ACC))]))[1], ]
#samp_meta[ACC == 'J012', UNI_ACC := ]

samp_meta[ACC == 'J421', UNI_ACC := 'SW787']
look_for_acc('SW787')

samp_meta[ACC == 'J422', UNI_ACC := 'SW788']
look_for_acc('SW788')

samp_meta[ACC == 'J423', UNI_ACC := 'SW789']
samp_meta[ACC == 'J423', STATE := 'Mississippi']
samp_meta[ACC == 'J423', LOCALITY := 'Multisite synthetic bred in MS; possibly used to generate SW790']
look_for_acc('SW789')

samp_meta[ACC == 'J424', UNI_ACC := 'SW790']
look_for_acc('SW790')

samp_meta[ACC == 'J425', UNI_ACC := 'SW793']
# Possibly PI659345
look_for_acc('SW793')

samp_meta[ACC == 'J426', UNI_ACC := 'SW795']
look_for_acc('SW795')

samp_meta[ACC == 'J427', UNI_ACC := 'SW796']
look_for_acc('SW796')

samp_meta[ACC == 'J428', UNI_ACC := 'SW797']
look_for_acc('SW797')

samp_meta[ACC == 'J429', UNI_ACC := 'SW798']
look_for_acc('SW798')

samp_meta[ACC == 'J430', UNI_ACC := 'SW799']
look_for_acc('SW799')

samp_meta[ACC == 'J431', UNI_ACC := 'SW802']
look_for_acc('SW802')

samp_meta[ACC == 'J432', UNI_ACC := 'SW803']
look_for_acc('SW803')

samp_meta[ACC == 'J433', UNI_ACC := 'SW805']
look_for_acc('SW805')

samp_meta[ACC == sort(unique(samp_meta$ACC[which(is.na(samp_meta$UNI_ACC))]))[1], ]
#samp_meta[ACC == 'J012', UNI_ACC := ]

samp_meta[ACC == 'J434', UNI_ACC := 'SW806']
# Might be PI 674680
look_for_acc('SW806')

samp_meta[ACC == 'J435', UNI_ACC := 'SW808']
look_for_acc('SW808')

samp_meta[ACC == 'J436', UNI_ACC := 'SW809']
look_for_acc('SW809')

samp_meta[ACC == 'J446', UNI_ACC := 'SWG32']
look_for_acc('SWG32')

samp_meta[ACC == 'J447', UNI_ACC := 'SWG39']
look_for_acc('SWG39')

samp_meta[ACC == 'J448', UNI_ACC := 'WS4U']
look_for_acc('WS4U')

samp_meta[ACC == 'J449', UNI_ACC := 'WS98_SB']
look_for_acc('WS98')

samp_meta[ACC == 'J450', UNI_ACC := 'ECS_1']
look_for_acc('ECS-1')

samp_meta[ACC == 'J451', UNI_ACC := 'ECS_2']
look_for_acc('ECS-2')

samp_meta[ACC == 'J454', UNI_ACC := 'ECS_12']
look_for_acc('ECS-12')

samp_meta[ACC == 'J455', UNI_ACC := 'ECS_6']
look_for_acc('ECS-6')

samp_meta[ACC == 'J456', UNI_ACC := 'GRIF_17004']
look_for_acc('2009FL')

samp_meta[ACC == 'J458', UNI_ACC := 'GRIF_17007']
look_for_acc('2009FL-045')

samp_meta[ACC == 'J460', UNI_ACC := 'GRIF_17504']
look_for_acc('17504')

samp_meta[ACC == 'J461', UNI_ACC := 'GRIF_17505']
look_for_acc('17505')

samp_meta[ACC == 'J462', UNI_ACC := 'GRIF_17511']
look_for_acc('17511')

samp_meta[ACC == 'J463', UNI_ACC := 'GRIF_17512']
look_for_acc('17512')

samp_meta[ACC == sort(unique(samp_meta$ACC[which(is.na(samp_meta$UNI_ACC))]))[1], ]
#samp_meta[ACC == 'J012', UNI_ACC := ]

samp_meta[ACC == 'J464', UNI_ACC := 'GRIF_17517']
look_for_acc('17517')

samp_meta[ACC == 'J465', UNI_ACC := 'GRIF_17521']
look_for_acc('17521')

samp_meta[ACC == 'J466', UNI_ACC := 'GRIF_17524']
look_for_acc('17524')

samp_meta[ACC == 'J467', UNI_ACC := 'GRIF_17537']
look_for_acc('17537')

samp_meta[ACC == 'J468', UNI_ACC := 'GRIF_17544']
look_for_acc('17544')

samp_meta[ACC == 'J469', UNI_ACC := 'GRIF_17005']
look_for_acc('17005')

samp_meta[ACC == 'J470', UNI_ACC := 'SRP']
look_for_acc('SRP')

samp_meta[ACC == 'J471', UNI_ACC := 'CHP']
look_for_acc('CHP')

samp_meta[ACC == 'J474', UNI_ACC := 'DWP']
look_for_acc('DWP')

samp_meta[ACC == 'J475', UNI_ACC := 'KYP']
look_for_acc('KYP')

samp_meta[ACC == 'J476', UNI_ACC := 'RRP']
look_for_acc('RRP')

samp_meta[ACC == 'J477', UNI_ACC := 'LCN']
look_for_acc('LCN')

samp_meta[ACC == 'J478', UNI_ACC := 'RTH']
look_for_acc('RTH')

samp_meta[ACC == 'J479', UNI_ACC := 'KGP']
look_for_acc('KGP')

samp_meta[ACC == 'J481', UNI_ACC := 'MBP']
look_for_acc('MBP')

samp_meta[ACC == 'J482', UNI_ACC := 'FLD']
look_for_acc('FLD')

samp_meta[ACC == sort(unique(samp_meta$ACC[which(is.na(samp_meta$UNI_ACC))]))[1], ]
#samp_meta[ACC == 'J012', UNI_ACC := ]

samp_meta[ACC == 'J483', UNI_ACC := 'ACR']
look_for_acc('ACR')

samp_meta[ACC == 'J484', UNI_ACC := 'PFL']
look_for_acc('PFL')

samp_meta[ACC == 'J485', UNI_ACC := 'SRG']
look_for_acc('SRG')

samp_meta[ACC == 'J486', UNI_ACC := 'HEF']
look_for_acc('HEF')

samp_meta[ACC == 'J488', UNI_ACC := 'DGP']
look_for_acc('DGP')

samp_meta[ACC == 'J489', UNI_ACC := 'OSP']
look_for_acc('OSP')

samp_meta[ACC == 'J490', UNI_ACC := 'WKT']
look_for_acc('WKT')

samp_meta[ACC == 'J491', UNI_ACC := 'TLL']
look_for_acc('TLL')

samp_meta[ACC == 'J494', UNI_ACC := 'TLD']
look_for_acc('TLD')

samp_meta[ACC == 'J495', UNI_ACC := 'EP']
look_for_acc('EP')

samp_meta[ACC == 'J496', UNI_ACC := 'GISP']
look_for_acc('GISP')

samp_meta[ACC == 'J497', UNI_ACC := 'HWY1']
look_for_acc('HWY1')

samp_meta[ACC == 'J498', UNI_ACC := 'BBM']
look_for_acc('BBM')

samp_meta[ACC == 'J499', UNI_ACC := 'FSP']
look_for_acc('FSP')

samp_meta[ACC == 'J500', UNI_ACC := 'BSP']
look_for_acc('BSP')

samp_meta[ACC == sort(unique(samp_meta$ACC[which(is.na(samp_meta$UNI_ACC))]))[1], ]
#samp_meta[ACC == 'J012', UNI_ACC := ]

samp_meta[ACC == 'J501', UNI_ACC := 'MP']
look_for_acc('MP')

samp_meta[ACC == 'J502', UNI_ACC := 'SSP']
look_for_acc('SSP')

samp_meta[ACC == 'J503', UNI_ACC := 'PBJSP']
look_for_acc('PBJSP')

samp_meta[ACC == 'J504', UNI_ACC := 'MC']
look_for_acc('MC')

samp_meta[ACC == 'J505', UNI_ACC := 'LBSP']
look_for_acc('LBSP')

samp_meta[ACC == 'J506', UNI_ACC := 'HP']
look_for_acc('HP')

samp_meta[ACC == 'J507', UNI_ACC := 'TSP']
look_for_acc('TSP')

samp_meta[ACC == 'J511', UNI_ACC := 'DUN']
look_for_acc('DUN')

samp_meta[ACC == 'J512', UNI_ACC := 'CLM']
look_for_acc('CLM')

samp_meta[ACC == 'J513', UNI_ACC := 'HEL']
look_for_acc('HEL')

samp_meta[ACC == 'J514', UNI_ACC := 'LRSP']
look_for_acc('JP 424')

samp_meta[ACC == 'J515', UNI_ACC := 'GRIF_16845']
look_for_acc('16845')

samp_meta[ACC == 'J516', UNI_ACC := 'GRIF_16853']
look_for_acc('16853')

samp_meta[ACC == 'J517', UNI_ACC := 'GRIF_16854']
look_for_acc('16854')

samp_meta[ACC == 'J518', UNI_ACC := 'GRIF_16884']
look_for_acc('16884')

samp_meta[ACC == sort(unique(samp_meta$ACC[which(is.na(samp_meta$UNI_ACC))]))[1], ]
#samp_meta[ACC == 'J012', UNI_ACC := ]

samp_meta[ACC == 'J520', UNI_ACC := 'PI_674658']
look_for_acc('674658')

samp_meta[ACC == 'J521', UNI_ACC := 'PI_674659']
look_for_acc('674659')

samp_meta[ACC == 'J522', UNI_ACC := 'PI_674660']
look_for_acc('674660')

samp_meta[ACC == 'J523', UNI_ACC := 'PI_674661']
look_for_acc('674661')

samp_meta[ACC == 'J525', UNI_ACC := 'PI_674663']
look_for_acc('674663')

samp_meta[ACC == 'J526', UNI_ACC := 'PI_674664']
look_for_acc('674664')

samp_meta[ACC == 'J527', UNI_ACC := 'PI_674665']
look_for_acc('674665')

samp_meta[ACC == 'J528', UNI_ACC := 'PI_674666']
look_for_acc('674666')

samp_meta[ACC == 'J529', UNI_ACC := 'PI_674667']
look_for_acc('674667')

samp_meta[ACC == 'J530', UNI_ACC := 'PI_674668']
look_for_acc('674668')

samp_meta[ACC == 'J531', UNI_ACC := 'PI_674669']
look_for_acc('674669')

samp_meta[ACC == 'J532', UNI_ACC := 'PI_674670']
look_for_acc('674670')

samp_meta[ACC == 'J533', UNI_ACC := 'PI_674671']
look_for_acc('674671')

samp_meta[ACC == 'J534', UNI_ACC := 'PI_674672']
look_for_acc('674672')

samp_meta[ACC == 'J535', UNI_ACC := 'Pederles']
look_for_acc('Pederles')

samp_meta[ACC == 'J536', UNI_ACC := 'PI_674674']
look_for_acc('674674')

samp_meta[ACC == sort(unique(samp_meta$ACC[which(is.na(samp_meta$UNI_ACC))]))[1], ]
#samp_meta[ACC == 'J012', UNI_ACC := ]

samp_meta[ACC == 'J537', UNI_ACC := 'PI_674676']
samp_meta[ACC == 'J538', UNI_ACC := 'PI_674676']
# based on locality notes, these are the same accession, though J537
#  *might* be PI 674675, since that also has the same locality notes
#  in GRIN and isn't in the list of accessions, that I can tell

samp_meta[ACC == 'J539', UNI_ACC := 'PI_674677']
look_for_acc('674677')

samp_meta[ACC == 'J540', UNI_ACC := 'PI_674678']
look_for_acc('674678')

samp_meta[ACC == 'J541', UNI_ACC := 'PI_674679']
look_for_acc('674679')

samp_meta[ACC == 'J542', UNI_ACC := 'PI_674680']
look_for_acc('674680')

samp_meta[ACC == 'J543', UNI_ACC := 'PI_674681']
look_for_acc('674681')

samp_meta[ACC == 'J544', UNI_ACC := 'PI_674682']
look_for_acc('674682')

samp_meta[ACC == 'J545', UNI_ACC := 'PI_674683']
look_for_acc('674683')

samp_meta[ACC == 'J546', UNI_ACC := 'PI_674684']
look_for_acc('674684')

samp_meta[ACC == 'J548', UNI_ACC := 'PI_674686']
look_for_acc('Delos')

samp_meta[ACC == 'J549', UNI_ACC := 'PI_674687']
look_for_acc('674687')

samp_meta[ACC == 'J554', UNI_ACC := 'PI_674692']
look_for_acc('674692')

samp_meta[ACC == 'J557', UNI_ACC := 'PI_674695']
look_for_acc('674695')

samp_meta[ACC == 'J558', UNI_ACC := 'PI_674696']
look_for_acc('674696')

samp_meta[ACC == sort(unique(samp_meta$ACC[which(is.na(samp_meta$UNI_ACC))]))[1], ]
#samp_meta[ACC == 'J012', UNI_ACC := ]

samp_meta[ACC == 'J559', UNI_ACC := 'PI_674697']
look_for_acc('674697')

samp_meta[ACC == 'J560', UNI_ACC := 'PI_674698']
look_for_acc('674698')

samp_meta[ACC == 'J561', UNI_ACC := 'PI_674699']
look_for_acc('674699')

samp_meta[ACC == 'J563', UNI_ACC := 'PI_674701']
look_for_acc('674701')

samp_meta[ACC == 'J564', UNI_ACC := 'PI_674702']
look_for_acc('674702')

samp_meta[ACC == 'J565', UNI_ACC := 'PI_674703']
look_for_acc('674703')

samp_meta[ACC == 'J566', UNI_ACC := 'PI_674704']
look_for_acc('674704')

samp_meta[ACC == 'J567', UNI_ACC := 'PI_674705']
look_for_acc('674705')

samp_meta[ACC == 'J568', UNI_ACC := 'PI_674706']
look_for_acc('674706')

samp_meta[ACC == 'J570', UNI_ACC := 'PI_674708']
look_for_acc('674708')

samp_meta[ACC == sort(unique(samp_meta$ACC[which(is.na(samp_meta$UNI_ACC))]))[1], ]
#samp_meta[ACC == 'J012', UNI_ACC := ]

samp_meta[ACC == 'J572', UNI_ACC := 'PI_674710']
look_for_acc('674710')

samp_meta[ACC == 'J574', UNI_ACC := 'PI_674712']
look_for_acc('674712')

samp_meta[ACC == 'J575', UNI_ACC := 'PI_674713']
look_for_acc('674713')

samp_meta[ACC == 'J576', UNI_ACC := 'GRIF_17453']
look_for_acc('17453')

samp_meta[ACC == 'J577', UNI_ACC := 'GRIF_17454']
look_for_acc('17454')

samp_meta[ACC == 'J578', UNI_ACC := 'GRIF_17458']
look_for_acc('17458')

samp_meta[ACC == 'J579', UNI_ACC := 'GRIF_17467']
look_for_acc('17467')

samp_meta[ACC == 'J580', UNI_ACC := 'FERMI']
look_for_acc('Fermi')

samp_meta[ACC == 'J581', UNI_ACC := 'JWP']
look_for_acc('JWP')

samp_meta[ACC == 'J583', UNI_ACC := 'AMGRD']
look_for_acc('Arcelor')

samp_meta[ACC == 'J584', UNI_ACC := 'HISP']
look_for_acc('HISP')

samp_meta[ACC == 'J585', UNI_ACC := 'PICP']
look_for_acc('PICP')

samp_meta[ACC == 'J586', UNI_ACC := 'HBSP']
look_for_acc('HBSP')

samp_meta[ACC == 'J587', UNI_ACC := 'HSF']
look_for_acc('HSF')

samp_meta[ACC == 'J588', UNI_ACC := 'NCBG_297']
look_for_acc('VA1')

samp_meta[ACC == sort(unique(samp_meta$ACC[which(is.na(samp_meta$UNI_ACC))]))[1], ]
#samp_meta[ACC == 'J012', UNI_ACC := ]

samp_meta[ACC == 'J589', UNI_ACC := 'NCBG_287']
look_for_acc('VA3')

samp_meta[ACC == 'J590', UNI_ACC := 'NCBG_541']
look_for_acc('VA4')

samp_meta[ACC == 'J591', UNI_ACC := 'NCBG_542']
look_for_acc('542')

samp_meta[ACC == 'J592', UNI_ACC := 'NY1']
look_for_acc('NY1')

samp_meta[ACC == 'J593', UNI_ACC := 'NCBG_276']
look_for_acc('NC1')

samp_meta[ACC == 'J594', UNI_ACC := 'NCBG_550']
look_for_acc('NC2')

samp_meta[ACC == 'J595', UNI_ACC := 'NC3']
look_for_acc('NC3')

samp_meta[ACC == 'J596', UNI_ACC := 'NCBG_529']
look_for_acc('MD1')

samp_meta[ACC == 'J597', UNI_ACC := 'NCBG_530']
look_for_acc('MD2')

samp_meta[ACC == 'J598', UNI_ACC := 'NCBG_490']
look_for_acc('MD3')

samp_meta[ACC == 'J599', UNI_ACC := 'NCBG_295']
look_for_acc('MD4')

samp_meta[ACC == 'J600', UNI_ACC := 'NCBG_540']
look_for_acc('MD5')

samp_meta[ACC == 'J601', UNI_ACC := 'NCBG_515']
look_for_acc('515')

samp_meta[ACC == 'J602', UNI_ACC := 'HBSC']
look_for_acc('Hobcaw')

samp_meta[ACC == sort(unique(samp_meta$ACC[which(is.na(samp_meta$UNI_ACC))]))[1], ]
#samp_meta[ACC == 'J012', UNI_ACC := ]

samp_meta[ACC == 'J603', UNI_ACC := 'FMSC']
look_for_acc('Francis')

samp_meta[ACC == 'J604', UNI_ACC := 'Patuxent']
look_for_acc('MD7')

samp_meta[ACC == 'J605', UNI_ACC := 'CT1']
look_for_acc('CT1')

samp_meta[ACC == 'J606', UNI_ACC := 'MA2']
look_for_acc('MA2')

samp_meta[ACC == 'J607', UNI_ACC := 'RI1']
look_for_acc('RI1')

samp_meta[ACC == 'J608', UNI_ACC := 'NWR']
look_for_acc('NWR')

samp_meta[ACC == 'J609', UNI_ACC := 'FMNF']
look_for_acc('SITE')

samp_meta[ACC == 'J610', UNI_ACC := 'KERR']
look_for_acc('JP77')

samp_meta[ACC == 'J612', UNI_ACC := 'PI_677185']
look_for_acc('677185')

samp_meta[ACC == 'J613', UNI_ACC := 'PI_476815']
look_for_acc('476815')

samp_meta[ACC == 'J614', UNI_ACC := 'PI_421136']
look_for_acc('421136')

samp_meta[ACC == 'J615', UNI_ACC := 'PI_561721']
look_for_acc('561721')
samp_meta[look_for_acc('561721'), TAXON := 'Panicum amarum']

samp_meta[ACC == 'J616', UNI_ACC := 'PI_677169']
look_for_acc('677169')

samp_meta[ACC == 'J617', UNI_ACC := 'PI_677149']
look_for_acc('677149')

samp_meta[ACC == 'J618', UNI_ACC := 'PI_677183']
look_for_acc('677183')

samp_meta[ACC == 'J619', UNI_ACC := 'PI_677181']
look_for_acc('677181')

samp_meta[ACC == sort(unique(samp_meta$ACC[which(is.na(samp_meta$UNI_ACC))]))[1], ]
#samp_meta[ACC == 'J012', UNI_ACC := ]

samp_meta[ACC == 'J620', UNI_ACC := 'NCBG_307']
look_for_acc('VA2')

samp_meta[ACC == 'J621', UNI_ACC := 'NCBG_594']
look_for_acc('MD8')

samp_meta[ACC == 'J623', UNI_ACC := 'TNFV']
look_for_acc('TNFV')

samp_meta[ACC == 'J625', UNI_ACC := 'PI_677149']
# based on locality info, looks like the same as J617
look_for_acc('Bluff Point')

samp_meta[ACC == 'J626', UNI_ACC := 'PI_677150']
# matched with GRIN PI number based on Locality info
look_for_acc('Lords Hill')

samp_meta[ACC == 'J627', UNI_ACC := 'PI_677153']
# matched with GRIN PI number based on Locality info
look_for_acc('Norwich')

samp_meta[ACC == 'J628', UNI_ACC := 'PI_677156']
# matched with GRIN PI number based on Locality info

samp_meta[ACC == 'J634', UNI_ACC := 'PI_677164']
# matched with GRIN PI number based on Locality info

samp_meta[ACC == 'J635', UNI_ACC := 'PI_677165']
# matched with GRIN PI number based on Locality info

samp_meta[ACC == 'J636', UNI_ACC := 'PI_677167']
# matched with GRIN PI number based on Locality info

samp_meta[ACC == 'J638', UNI_ACC := 'PI_677173']
# matched with GRIN PI number based on Locality info

samp_meta[ACC == 'J639', UNI_ACC := 'PI_677174']
# matched with GRIN PI number based on Locality info

samp_meta[ACC == 'J640', UNI_ACC := 'PI_677175']
# matched with GRIN PI number based on Locality info

samp_meta[ACC == 'J641', UNI_ACC := 'PI_677176']
# matched with GRIN PI number based on Locality info

samp_meta[ACC == 'J643', UNI_ACC := 'PI_677180']
# matched with GRIN PI number based on Locality info

samp_meta[ACC == sort(unique(samp_meta$ACC[which(is.na(samp_meta$UNI_ACC))]))[1], ]
#samp_meta[ACC == 'J012', UNI_ACC := ]

samp_meta[ACC == 'J644', UNI_ACC := 'PI_677182']
# matched with GRIN PI number based on Locality info

samp_meta[ACC == 'J645', UNI_ACC := 'PI_677183']
# matched with GRIN PI number based on Locality info

samp_meta[ACC == 'J646', UNI_ACC := 'PI_677184']
# matched with GRIN PI number based on Locality info

samp_meta[ACC == 'J647', UNI_ACC := 'PI_677185']
# matched with GRIN PI number based on Locality info

samp_meta[ACC == 'J649', UNI_ACC := 'LINP']
look_for_acc('Native Plants')

samp_meta[ACC == 'J651', UNI_ACC := 'COLONY']
look_for_acc('Colony')

samp_meta[ACC == 'J652', UNI_ACC := 'BOMASTER']
look_for_acc('BoMaster')

#samp_meta[ACC == 'J652', UNI_ACC := []
samp_meta[look_for_acc('Stuart'), UNI_ACC := 'STUART']
look_for_acc('Stuart')

samp_meta[ACC == 'J656', UNI_ACC := 'DOL']
look_for_acc('DOL')

samp_meta[ACC == 'J657', UNI_ACC := 'IGN']
look_for_acc('IGN')

samp_meta[ACC == 'J658', UNI_ACC := 'DAN']
look_for_acc('DAN')

samp_meta[ACC == 'J659', UNI_ACC := 'LBJNG']
look_for_acc('LBJ')

samp_meta[ACC == 'J661', UNI_ACC := 'VA6']
look_for_acc('VA6')

samp_meta[ACC == 'J662', UNI_ACC := 'VA7']
look_for_acc('VA7')

samp_meta[ACC == 'J663', UNI_ACC := 'MD9']
look_for_acc('MD9')

samp_meta[ACC == 'J664', UNI_ACC := 'NJ1']
look_for_acc('NJ1')

samp_meta[ACC == sort(unique(samp_meta$ACC[which(is.na(samp_meta$UNI_ACC))]))[1], ]
#samp_meta[ACC == 'J012', UNI_ACC := ]

samp_meta[ACC == 'J665', UNI_ACC := 'NY2']
look_for_acc('NY2')

samp_meta[ACC == 'J666', UNI_ACC := 'NC4']
look_for_acc('NC4')

samp_meta[ACC == 'J667', UNI_ACC := 'NC5']
look_for_acc('NC5')

samp_meta[ACC == 'J668', UNI_ACC := 'NC6']
look_for_acc('NC6')

samp_meta[ACC == 'J669', UNI_ACC := 'NY3']
look_for_acc('NY3')

samp_meta[ACC == 'J670', UNI_ACC := 'NJ2']
look_for_acc('NJ2')

samp_meta[ACC == 'J671', UNI_ACC := 'MD10']
look_for_acc('MD10')

samp_meta[ACC == 'J672', UNI_ACC := 'DE2']
look_for_acc('DE2')

samp_meta[ACC == 'J673', UNI_ACC := 'DE3']
look_for_acc('DE3')

samp_meta[ACC == 'J674', UNI_ACC := 'DE1']
look_for_acc('DE1')

samp_meta[ACC == 'J675', UNI_ACC := 'VA10']
look_for_acc('VA10')

samp_meta[ACC == 'J676', UNI_ACC := 'VA9']
look_for_acc('VA9')

samp_meta[ACC == 'J677', UNI_ACC := 'NC7']
look_for_acc('NC7')

samp_meta[ACC == 'J678', UNI_ACC := 'NC8']
look_for_acc('NC8')

samp_meta[ACC == 'J679', UNI_ACC := 'NC9']
look_for_acc('NC9')

samp_meta[ACC == 'J680', UNI_ACC := 'NC10']
look_for_acc('NC10')

samp_meta[ACC == sort(unique(samp_meta$ACC[which(is.na(samp_meta$UNI_ACC))]))[1], ]
#samp_meta[ACC == 'J012', UNI_ACC := ]

samp_meta[ACC == 'J681', UNI_ACC := 'THUNDERCLOUD']
look_for_acc('Thunder')

samp_meta[ACC == 'J682', UNI_ACC := 'NJ4']
look_for_acc('NJ4')

samp_meta[ACC == 'J683', UNI_ACC := 'NJ3']
look_for_acc('NJ3')

samp_meta[ACC == 'J684', UNI_ACC := 'DE4']
look_for_acc('DE4')

samp_meta[ACC == 'J685', UNI_ACC := 'DE5']
look_for_acc('DE5')

samp_meta[which(is.na(samp_meta$UNI_ACC))[1]]

mx_na_inds <- intersect(which(is.na(samp_meta$UNI_ACC)), 
  grep('MX_', samp_meta$VCF_NAME))

samp_meta$UNI_ACC[mx_na_inds] <- samp_meta$VCF_NAME[mx_na_inds]

samp_meta[VCF_NAME == 'WBC3a', UNI_ACC := 'WBC']
look_for_acc('Bee Caves')

samp_meta[VCF_NAME == 'A4', UNI_ACC := 'ALAMO']
look_for_acc('A4')

#samp_meta[VCF_NAME == 'NFGA37_05', UNI_ACC := ]
samp_meta[look_for_acc('EG1102'), UNI_ACC := 'EG1102']
look_for_acc('EG1102')

samp_meta[VCF_NAME == 'NFGA09_05', UNI_ACC := 'KANLOW']

samp_meta[VCF_NAME == 'NFGA02_06', UNI_ACC := 'PI_414065']
look_for_acc('414065')

samp_meta[VCF_NAME == 'Hammonasset', UNI_ACC := 'HAMMONASSET']
look_for_acc('Hammonasset')

samp_meta[VCF_NAME == 'NL94', UNI_ACC := 'KANLOW_LIKE_NL94']

samp_meta[VCF_NAME == 'Pv346_1', UNI_ACC := 'AP13_VS16_CROSS']
samp_meta[VCF_NAME == 'Pv458_1', UNI_ACC := 'AP13_VS16_CROSS']

samp_meta[VCF_NAME == 'Dacotah_WO1', UNI_ACC := 'DACOTAH']

samp_meta[VCF_NAME == 'Pv304_1', UNI_ACC := 'AP13_VS16_CROSS']
samp_meta[VCF_NAME == 'Pv317_1', UNI_ACC := 'AP13_VS16_CROSS']

samp_meta[which(is.na(samp_meta$UNI_ACC))[1]]

samp_meta[VCF_NAME == '9001-3_BN389-69S', UNI_ACC := 'BN389']
samp_meta[VCF_NAME == '9020-5', UNI_ACC := 'CARTHAGE']

samp_meta[VCF_NAME == 'SW786_WO1', UNI_ACC := 'CONTAMINANT_1']

samp_meta[VCF_NAME == '2785b', UNI_ACC := 'KANLOW_X_SUMMER']

samp_meta[VCF_NAME == 'NL94_85_1', UNI_ACC := 'KANLOW_LIKE_NL94']
samp_meta[VCF_NAME == 'NL94_85_1_5', UNI_ACC := 'KANLOW_LIKE_NL94']

samp_meta[grep('SL93', samp_meta$VCF_NAME), UNI_ACC := 'ALAMO_LIKE_SL93']

samp_meta[PLANT_ID == 'AP13', UNI_ACC := 'AP13']

# look at patterns

table(table(samp_meta$UNI_ACC))
#   1   2   3   4   5   6   7   8  13  21  39 
# 213 104 132  20   6   4   2   1   2   1   1 
# ND samples from Arvid Boe have 39 and 21 samples;
# Alamo and Kanlow each have 13 samples

fwrite(samp_meta, out_file, na = 'NA')

quit(save = 'no')

#tmp_workspace <- '/home/f2p1/work/grabowsk/data/switchgrass/sg_8X_metadata_v1.1_tmp_workspace.Rdata'
#save.image(file = tmp_workspace)

