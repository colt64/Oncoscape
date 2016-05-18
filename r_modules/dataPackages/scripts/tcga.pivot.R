#Reference Table#------------------------------------------------------------------------------------
options(stringsAsFactors=FALSE)
library(RUnit)
library(R.utils)
library(stringr)
library(plyr)

## store dxyear
#  'initial_pathologic_dx_year' = list(name = "dxyear", data = "tcgaDate")


#--------------------------------------------------------------------------------
date.Calculation <- function(from_date, days){
	date <- format(as.Date(from_date, "%m/%d/%Y") + as.integer(days), "%m/%d/%Y")
	return(date)
}	



#----------------------     DOB functions Start Here      -----------------------
	BIRTH.unique.request <- function(study_name){
	  barcode
	  unique.dob <- unique(df$dob)
	  unique.gender <- unique(df$gender)
	  unique.ethnicity <- unique(df$ethnicity)
	  unique.race <- unique(df$race)
	  result = list(unique.dob=unique.dob, unique.gender=unique.gender, 
	  				unique.ethnicity=unique.ethnicity, unique.race=unique.race)
	  return(result)
	}

#----------------------   Diagnosis functions Start Here   ----------------------
	Diagnosis.unique.request <- function(study_name){
	 'bcr_patient_barcode' = list(name = "PatientID", data = "tcgaId"),
	 'tumor_tissue_site' = list(name = "disease", data ="upperCharacter"),
	 'tissue_source_site' = list(name = "tissueSourceSiteCode", data = "upperCharacter"),
	 'initial_pathologic_dx_year' = list(name = "dxyear", data = "tcgaDate")
   
	  unique.disease <- unique(df$disease)
	  unique.tissueSourceSiteCode <- unique(df$tissueSourceSiteCode)
	  result = list(unique.disease=unique.disease, 
	  				unique.tissueSourceSiteCode=unique.tissueSourceSiteCode)
	  return(result)
	}
#----------------------   Drug functions Start Here   ---------------------------
	Drug.unique.request <- function(study_name){
		tbl.drug 
			 'bcr_patient_barcode' = list(name = "PatientID", data = "tcgaId"),
			 'pharmaceutical_tx_started_days_to' = list(name = "drugStart", data = "character"),
			 'pharmaceutical_tx_ended_days_to' = list(name = "drugEnd", data = "character"),
			 'pharmaceutical_therapy_drug_name' = list(name = "agent", data = "upperCharacter"),
			 'pharmaceutical_therapy_type' = list(name = "therapyType", data = "upperCharacter"),
			 'therapy_regimen' = list(name = "intent", data = "upperCharacter"),
			 'prescribed_dose' = list(name = "dose", data = "upperCharacter"),
			 'total_dose' = list(name = "totalDose", data = "upperCharacter"),
			 'pharmaceutical_tx_dose_units' = list(name = "units", data = "upperCharacter"),
			 'pharmaceutical_tx_total_dose_units' = list(name = "totalDoseUnits", data = "upperCharacter"),
			 'route_of_administration' = list(name = "route", data = "upperCharacter"),
			 'pharma_adjuvant_cycles_count' = list(name = "cycle", data = "upperCharacter")
		
		tbl.omf 
			 'bcr_patient_barcode' = list(name = "PatientID", data = "tcgaId"),
			 'drug_name' = list(name = "agent", data = "upperCharacter"),
			 'days_to_drug_therapy_start' = list(name = "drugStart", data = "character"),
			 'malignancy_type' = list(name = "intent", data = "upperCharacter")
		
	    # reorganize three tbls 
	    tbl.f <- data.frame()
		tbl.f <- rbind.fill(tbl.f, tbl.drug)
		tbl.f <- rbind.fill(tbl.f, tbl.omf)

	    	data.Chemo <- merge(tbl.drug, tbl.pt, by = "PatientID", all.x = T)
		    data.Chemo$start <- rep(NA,nrow(data.Chemo))
		    data.Chemo$end <- rep(NA,nrow(data.Chemo))
		  	
		  	df <- data.Chemo
		  	unique.drugStart <- unique(df$drugStart)
		  	unique.drugEnd <- unique(df$drugEnd)
			unique.therapyType <- unique(df$therapyType)
			unique.intent <- unique(df$intent)
			unique.dose <- unique(df$dose)
			unique.units <- unique(df$units)
			unique.totalDose <- unique(df$totalDose)
			unique.totalDoseUnits <- unique(df$totalDoseUnits)
			unique.route <- unique(df$route)
			unique.cycle <- unique(df$cycle)
		  	result = list(unique.drugStart=unique.drugStart, 
						  unique.drugEnd=unique.drugEnd, 
		  				  unique.therapyType=unique.therapyType, 
		  				  unique.intent=unique.intent,
		  				  unique.dose=unique.dose,
		  				  unique.units=unique.units,
		  				  unique.totalDose=unique.totalDose,
		  				  unique.totalDoseUnits=unique.totalDoseUnits,
		  				  unique.route=unique.route,
		  				  unique.cycle=unique.cycle)
		  	print(study_name)
		  	return(result)
	    }
	}        
	#--------------------------------------------------------------------------------
	Drug.mapping.date <- function(df){
		df$start[which(is.na(df$drugStart))] <- NA
		df$end[which(is.na(df$drugEnd))] <- NA
		df[which(is.na(df$dxyear)), c("drugStart","drugEnd")] <- NA
	   
	    df$start <- format(as.Date(df$dxyear, "%m/%d/%Y") + as.integer(df$drugStart), "%m/%d/%Y")
	    df$end <- format(as.Date(df$dxyear, "%m/%d/%Y") + as.integer(df$drugEnd), "%m/%d/%Y")
			
		return(df)
	}	
	#--------------------------------------------------------------------------------
	Drug.mapping.filling <- function(df){
		requiredFields <- c("PatientID", "start", "end", "agent", "therapyType", "intent", "dose",
							"units", "totalDose", "totalDoseUnits", "route", "cycle")   				
		m <- matrix(nrow=nrow(df), ncol=length(which(!requiredFields %in% colnames(df))))
	    m <- as.data.frame(m)
	    colnames(m) <- requiredFields[(which(!(requiredFields) %in% colnames(df)))]
	    df <- cbind(df, m)     				
		return(df)
	}		
#----------------------   Radiation functions Start Here   ----------------------
	Rad.unique.request <- function(study_name){
		tbl.rad 
			 'bcr_patient_barcode' = list(name = "PatientID", data = "tcgaId"),
			 'radiation_therapy_started_days_to' = list(name = "radStart", data = "character"),
			 'radiation_therapy_ended_days_to' = list(name = "radEnd", data = "character"),
			 'radiation_therapy_type' = list(name = "radType", data = "upperCharacter"),
			 'radiation_type_other' = list(name = "radTypeOther", data = "upperCharacter"),
			 'therapy_regimen' = list(name = "intent", data = "upperCharacter"),
			 'radiation_therapy_site' = list(name = "target", data = "upperCharacter"),
			 'radiation_total_dose' = list(name = "totalDose", data = "upperCharacter"),
			 'radiation_adjuvant_units' = list(name = "totalDoseUnits", data = "upperCharacter"),
			 'radiation_adjuvant_fractions_total' = list(name = "numFractions", data = "upperCharacter")
		tbl.omf 
			 'bcr_patient_barcode' = list(name = "PatientID", data = "tcgaId"),
			 'radiation_tx_extent' = list(name = "target", data = "upperCharacter"),
			 'rad_tx_to_site_of_primary_tumor' = list(name = "targetAddition", data = "upperCharacter"),
			 'days_to_radiation_therapy_start' = list(name = "radStart", data = "character")

		tbl.f <- data.frame()
		tbl.f <- rbind.fill(tbl.f, tbl.rad)
		tbl.f <- rbind.fill(tbl.f, tbl.omf)

	    	data.Rad <- merge(tbl.f, tbl.pt, by = "PatientID", all.x = T)
		    data.Rad$start <- rep(NA,nrow(data.Rad))
		    data.Rad$end <- rep(NA,nrow(data.Rad))
		  	
		  	df <- data.Rad
		  	unique.radStart <- unique(df$radStart)
		  	unique.radEnd <- unique(df$radEnd)
			unique.radType <- unique(df$radType)
			unique.radTypeOther <- unique(df$radTypeOther)
			unique.intent <- unique(df$intent)
			unique.target <- unique(df$target)
			unique.targetAddition <- unique(df$targetAddition)
			unique.totalDose <- unique(df$totalDose)
			unique.totalDoseUnits <- unique(df$totalDoseUnits)
			unique.numFractions <- unique(df$numFractions)
		  	result = list(unique.radStart=unique.radStart, 
						  unique.radEnd=unique.radEnd, 
		  				  unique.radType=unique.radType,
		  				  unique.radTypeOther=unique.radTypeOther, 
		  				  unique.intent=unique.intent,
		  				  unique.target=unique.target,
		  				  unique.targetAddition=unique.targetAddition,
		  				  unique.totalDose=unique.totalDose,
		  				  unique.totalDoseUnits=unique.totalDoseUnits,
		  				  unique.numFractions=unique.numFractions)
		  	print(study_name)
		  	return(result)
	#--------------------------------------------------------------------------------	   
	    df$start <- format(as.Date(df$dxyear, "%m/%d/%Y") + as.integer(df$radStart), "%m/%d/%Y")
	    df$end <- format(as.Date(df$dxyear, "%m/%d/%Y") + as.integer(df$radEnd), "%m/%d/%Y")
	#--------------------------------------------------------------------------------

		tmpRadType <-rad_ref[match(df$radType,rad_ref$COMMON.RADTYPE),]$STANDARDIZED.RADTYPE
		tmpRadTypeOther <-rad_ref[match(df$radTypeOther,rad_ref$COMMON.OTHERRADTYPE),]$STANDARDIZED.OTHERRADTYPE
		tmpRadType[which(tmpRadType == "OTHER: SPECIFY IN NOTES")] <- tmpRadTypeOther[which(tmpRadType == "OTHER: SPECIFY IN NOTES")]

		df$radType <- tmpRadType
		return(df)
	}	
	#--------------------------------------------------------------------------------
		df$totalDoseUnits[grep("CGY", df$totalDose)] <- "CGY"
		df$totalDose[grep("CGY", df$totalDose)] <- str_extract(df$totalDose[grep("CGY", df$totalDose)],"[0-9,]+")
		
		df$totalDose[which(df$totalDoseUnits == "GY")] <- 
					as.integer(df$totalDose[which(df$totalDoseUnits == "GY")]) * 100

		df$totalDoseUnits[which(df$totalDoseUnits == "GY")] <- "CGY"
		return(df)
	}	
	#--------------------------------------------------------------------------------
	Rad.mapping.filling <- function(df){	
		requiredFields <- c("PatientID", "start", "end", "radType", "intent", "target",
							 "totalDose", "totalDoseUnits", "numFractions")   				
		m <- matrix(nrow=nrow(df), ncol=length(which(!requiredFields %in% colnames(df))))
	    m <- as.data.frame(m)
	    colnames(m) <- requiredFields[(which(!(requiredFields) %in% colnames(df)))]
	    df <- cbind(df, m)     				
		return(df)
	}			
#----------------------     Status functions Start Here      --------------------
	Status.unique.request <- function(study_name){
		tbl.pt 
			 'bcr_patient_barcode' = list(name = "PatientID", data = "tcgaId"),
			 'initial_pathologic_dx_year' = list(name = "dxyear", data = "tcgaDate"),
			 'vital_status' = list(name = "vital", data = "upperCharacter"),
			 'tumor_status' = list(name = "tumorStatus", data = "upperCharacter"),
			 'last_contact_days_to' = list(name = "lastContact", data = "upperCharacter"),
			 'death_days_to' = list(name = "deathDate", data = "upperCharacter")
		tbl.f1 
			 'bcr_patient_barcode' = list(name = "PatientID", data = "tcgaId"),
			 'vital_status' = list(name = "vital", data = "upperCharacter"),
			 'tumor_status' = list(name = "tumorStatus", data = "upperCharacter"),
			 'last_contact_days_to' = list(name = "lastContact", data = "upperCharacter"),
			 'death_days_to' = list(name = "deathDate", data = "upperCharacter")
		tbl.f2 
			 'bcr_patient_barcode' = list(name = "PatientID", data = "tcgaId"),
			 'vital_status' = list(name = "vital", data = "upperCharacter"),
			 'tumor_status' = list(name = "tumorStatus", data = "upperCharacter"),
			 'last_contact_days_to' = list(name = "lastContact", data = "upperCharacter"),
			 'death_days_to' = list(name = "deathDate", data = "upperCharacter")
		tbl.f3 
			 'bcr_patient_barcode' = list(name = "PatientID", data = "tcgaId"),
			 'vital_status' = list(name = "vital", data = "upperCharacter"),
			 'tumor_status' = list(name = "tumorStatus", data = "upperCharacter"),
			 'last_contact_days_to' = list(name = "lastContact", data = "upperCharacter"),
			 'death_days_to' = list(name = "deathDate", data = "upperCharacter")

		tbl.f <- tbl.pt[,-grep("dxyear", colnames(tbl.pt))]
		tbl.f <- rbind.fill(tbl.f, tbl.f1)
		tbl.f <- rbind.fill(tbl.f, tbl.f2)
		tbl.f <- rbind.fill(tbl.f, tbl.f3)
		df <- tbl.f
		unique.deathDate <- unique(df$deathDate)
		unique.lastContact <- unique(df$lastContact)
	  	unique.vital <- unique(df$vital)
	  	if("tumorStatus" %in% colnames(df)){
	  		unique.tumorStatus <- unique(df$tumorStatus)
	  		}else{
	  			unique.tumorStatus <- NA
	  		}
	  	
	  	result = list(unique.deathDate=unique.deathDate, unique.lastContact=unique.lastContact, 
	  				  unique.vital=unique.vital, unique.tumorStatus=unique.tumorStatus)
	  	return(result)
	}
	#--------------------------------------------------------------------------------
	Status.unique.aggregate <- function(res1, res2){
		res = list(unique.deathDate=unique(c(res1$unique.deathDate,res2$unique.deathDate)),
				   unique.lastContact=unique(c(res1$unique.lastContact,res2$unique.lastContact)),
				   unique.vital=unique(c(res1$unique.vital,res2$unique.vital)),
				   unique.tumorStatus=unique(c(res1$tumorStatus.race, res2$unique.tumorStatus)))
	    return(res)
	}
	#--------------------------------------------------------------------------------
	Status.mapping.date.Check <- function(df){
		
		if(length(which(as.numeric(df$lastContact) > as.numeric(df$deathDate)))) {
			lastContactGreaterThanDeath  = paste(df[which(df$lastContact > df$deathDate),]$PatientID)
			warning("last contact occured after death: ", lastContactGreaterThanDeath)
		}
        df[which(!(is.na(df$lastContact))),]$date <- df[which(!(is.na(df$lastContact))),]$lastContact
        df[which(!(is.na(df$deathDate))),]$date <- df[which(!(is.na(df$deathDate))),]$deathDate

		return(df)
	}	
	#--------------------------------------------------------------------------------
	Status.mapping.filling <- function(df){	
		requiredFields <- c("PatientID", "date", "vital", "tumorStatus")   				
		m <- matrix(nrow=nrow(df), ncol=length(which(!requiredFields %in% colnames(df))))
	    m <- as.data.frame(m)
	    colnames(m) <- requiredFields[(which(!(requiredFields) %in% colnames(df)))]
	    df <- cbind(df, m)     				
		return(df)
	}			
#----------------------   PROGRESSION functions Start Here   ----------------------
	Progression.unique.request <- function(study_name){
		tbl.f1 
		 'bcr_patient_barcode' = list(name = "PatientID", data = "tcgaId"),
		 'new_tumor_event_dx_days_to' = list(name = "newTumorDate", data = "upperCharacter"),
		 'new_neoplasm_event_type' = list(name = "newTumor", data = "upperCharacter"),
		 'new_tumor_event_type' = list(name = "newTumor", data = "upperCharacter")
		tbl.f2 
		 'bcr_patient_barcode' = list(name = "PatientID", data = "tcgaId"),
		 'new_tumor_event_dx_days_to' = list(name = "newTumorDate", data = "upperCharacter"),
		 'new_neoplasm_event_type' = list(name = "newTumor", data = "upperCharacter"),
		 'new_tumor_event_type' = list(name = "newTumor", data = "upperCharacter")
		tbl.nte 
		 'bcr_patient_barcode' = list(name = "PatientID", data = "tcgaId"),
		 'new_tumor_event_dx_days_to' = list(name = "newTumorDate", data = "upperCharacter"),
		 'new_neoplasm_event_type' = list(name = "newTumor", data = "upperCharacter"),
		 'new_tumor_event_type' = list(name = "newTumor", data = "upperCharacter")
		tbl.nte_f1 
		 'bcr_patient_barcode' = list(name = "PatientID", data = "tcgaId"),
		 'new_tumor_event_dx_days_to' = list(name = "newTumorDate", data = "upperCharacter"),
		 'new_neoplasm_event_type' = list(name = "newTumor", data = "upperCharacter"),
		 'new_tumor_event_type' = list(name = "newTumor", data = "upperCharacter")
		tbl.f <- data.frame()
		if(exists("tbl.nte")) tbl.f <- rbind.fill(tbl.f, tbl.nte)
		if(exists("tbl.f1")) tbl.f <- rbind.fill(tbl.f, tbl.f1)
		if(exists("tbl.f2")) tbl.f <- rbind.fill(tbl.f, tbl.f2)
		if(exists("tbl.nte_f1")) tbl.f <- rbind.fill(tbl.f, tbl.nte_f1)
			df <- merge(tbl.f, tbl.pt)
			unique.newTumor <- unique(df$newTumor)
			unique.newTumorDate <- unique(df$newTumorDate)
		   	result = list(unique.newTumor=unique.newTumor, unique.newTumorDate=unique.newTumorDate)
		  	return(result)
		}
	#--------------------------------------------------------------------------------
	Progression.mapping.newTumorNaRM <- function(df){
		rmList <- apply(df, 1, function(x){
					pt = getElement(x, "PatientID")
					newTumorType = getElement(x, "newTumor")
					newTumorDateVal = getElement(x, "newTumorDate")
					tmp <- subset(df, PatientID == pt & newTumorDate == newTumorDateVal)
					if(nrow(tmp) > 1 && any(is.na(tmp$newTumor))){
						return(which(df$PatientID == pt & df$newTumorDate == newTumorDateVal & is.na(df$newTumor)))
					}
				})
		if(is.null(rmList)){
			return(df)
		}else{
			#print(unlist(rmList))
			df <- df[-(unlist(rmList)),]
		}
	}
#----------------------   Encounter functions Start Here   ----------------------
  # brca, hnsc, prad DO NOT HAVE ENCOUNTER RECORDS!
  Encounter.unique.request <- function(study_name){   
	    #(tbl.pt 'encType','karnofsky_score','ECOG only in gbm,lgg,luad,lusc)
	    tbl.pt 
		 'bcr_patient_barcode' = list(name = "PatientID", data = "tcgaId"),
		 'performance_status_timing' = list(name = "encType", data = "upperCharacter"),
		 'karnofsky_score'= list(name = "KPS", data = "upperCharacter"),
		 'ecog_score' = list(name = "ECOG", data = "upperCharacter"),
		 #coad/read only
		 'height_cm_at_diagnosis' = list(name = "height", data = "upperCharacter"),
		 'weight_kg_at_diagnosis' = list(name = "weight", data = "upperCharacter"),
		 #lung only
		 'fev1_fvc_ratio_prebroncholiator'= list(name = "prefev1.ratio", data = "upperCharacter"),
		 'fev1_percent_ref_prebroncholiator'= list(name = "prefev1.percent", data = "upperCharacter"),
		 'fev1_fvc_ratio_postbroncholiator'= list(name = "postfev1.ratio", data = "upperCharacter"),
		 'fev1_percent_ref_postbroncholiator'= list(name = "postfev1.percent", data = "upperCharacter"),
		 'carbon_monoxide_diffusion_dlco'= list(name = "carbon.monoxide.diffusion", data = "upperCharacter")

	    #(tbl.f1'encType','karnofsky_score','ECOG only in gbm,lgg,luad,lusc)
	    tbl.f1 
		   list('bcr_patient_barcode' = list(name = "PatientID", data = "tcgaId"),
				'performance_status_timing' = list(name = "encType", data = "upperCharacter"),
				'karnofsky_score'= list(name = "KPS", data = "upperCharacter"),
				'ecog_score' = list(name = "ECOG", data = "upperCharacter")

	    data.Encounter <- rbind.fill(tbl.pt, tbl.f1)
	    #colnames(data.Encounter)
	    
	    df <- data.Encounter
	    unique.encType<- unique(df$encType)
	    unique.KPS <- unique(df$KPS)
	    unique.ECOG <- unique(df$ECOG)
	    #coad/read only
	    unique.height <- unique(df$height)
	    unique.weight <- unique(df$weight)
	    #lung only
	    unique.prefev1.ratio <- unique(df$prefev1.ratio)
	    unique.prefev1.percent <- unique(df$prefev1.percent)
	    unique.postfev1.ratio<- unique(df$postfev1.ratio)
	    unique.postfev1.percent <- unique(df$postfev1.percent)
	    unique.carbon.monoxide.diffusion<- unique(df$carbon.monoxide.diffusion)
	    
	    
	    result = list(unique.encType=unique.encType, 
	                  unique.KPS=unique.KPS,
	                  unique.ECOG=unique.ECOG,
	                  unique.height=unique.height,
	                  unique.weight=unique.weight,
	                  unique.prefev1.ratio=unique.prefev1.ratio,
	                  unique.prefev1.percent=unique.prefev1.percent,
	                  unique.postfev1.ratio=unique.postfev1.ratio,
	                  unique.postfev1.percent=unique.postfev1.percent,
	                  unique.carbon.monoxide.diffusion=unique.carbon.monoxide.diffusion)
	    print(study_name)
	    return(result)
  }
#----------------------   Procedure functions Start Here   ----------------------
  	Procedure.unique.request <- function(study_name){
	    tbl.nte 
		 'bcr_patient_barcode' = list(name = "PatientID", data = "tcgaId"),
		 'new_tumor_event_surgery_days_to_loco' = list(name = "date_locoregional", data = "upperCharacter"), #(only in lgg,luad,lusc)
		 'new_tumor_event_surgery_days_to_met'= list(name = "date_metastatic", data = "upperCharacter"), #(only in lgg,luad,lusc)
		 #'new_tumor_event_surgery' = list(name = "new_tumor_event_surgery", data = "upperCharacter"), #(in brca,hnsc but not being collected...) #YES/NO
		 'days_to_new_tumor_event_additional_surgery_procedure'  = list(name = "date", data = "upperCharacter"), #(only in gbm,coad,read)
		 'new_neoplasm_event_type'  = list(name = "site", data = "upperCharacter"), #(only in gbm, coad, read)
		 'new_tumor_event_type'  = list(name = "site", data = "upperCharacter") #(only in hnsc, prad, luad, lusc)
		 #'new_tumor_event_additional_surgery_procedure'  = list(name = "new_tumor_event_additional_surgery_procedure", data = "upperCharacter") #(gbm,coad,read but not being collected...) YES/NO
	    tbl.omf
			 'bcr_patient_barcode' = list(name = "PatientID", data = "tcgaId"),
			 'days_to_surgical_resection' = list(name = "date", data = "upperCharacter"), #(gbm,lgg,hnsc,brca,prad,luad,lusc,coad,read)
			 'other_malignancy_laterality' = list(name = "side", data = "upperCharacter"), #(brca)
			 'surgery_type' = list(name = "surgery_name", data = "upperCharacter") #(gbm,lgg,hnsc,brca,pProcedure,lusc,luad,coad,read) 
	    tbl.pt 
			   'bcr_patient_barcode' = list(name = "PatientID", data = "tcgaId"),
			   'initial_pathologic_dx_year' = list(name = "dxyear", data = "tcgaDate"),
			   'laterality'  = list(name = "side", data = "upperCharacter"), #(only in lgg, hnsc, prad)
			   'tumor_site' = list(name = "site", data = "upperCharacter"),  #(only in lgg)
			   'supratentorial_localization'= list(name = "site", data = "upperCharacter"), #(only in lgg)
			   'surgical_procedure_first'= list(name = "surgery_name", data = "upperCharacter"), #only in brca
			   'first_surgical_procedure_other'= list(name = "surgery_name", data = "upperCharacter") #only in brca
	    tbl.f1 
			   'bcr_patient_barcode' = list(name = "PatientID", data = "tcgaId"),
			   'new_tumor_event_surgery_days_to_loco' = list(name = "date_locoregional", data = "upperCharacter"), #(only in lgg,hnsc,luad,lusc)
			   'new_tumor_event_surgery_days_to_met'= list(name = "date_metastatic", data = "upperCharacter") #(only in lgg,hnsc,luad,lusc)
			   #'new_tumor_event_surgery' = list(name = "new_tumor_event_surgery", data = "upperCharacter") #(In lgg,luad,lusc but not being collected...)
	 
	      tbl.nte_f1 
			   'bcr_patient_barcode' = list(name = "PatientID", data = "tcgaId"),
			   #'new_tumor_event_surgery' = list(name = "new_tumor_event_surgery", data = "upperCharacter"), #(used to build hnsc tables but is also a column in brca that is not being collected)
			   'days_to_new_tumor_event_additional_surgery_procedure'  = list(name = "date", data = "upperCharacter"), #(only in gbm,hnsc,coad,read)
			   'new_neoplasm_event_type'  = list(name = "site", data = "upperCharacter"), #(only in gbm, coad, read)
			   'new_tumor_event_type'  = list(name = "site", data = "upperCharacter") #(only in hnsc, brca)
			   #'new_tumor_event_additional_surgery_procedure'  = list(name = "new_tumor_event_additional_surgery_procedure", data = "upperCharacter") #(hnsc)
	
		data.Procedure <- data.frame()
			data.Procedure <- rbind.fill(tbl.pt[, -match("dxyear", colnames(tbl.pt))], data.Procedure)
			data.Procedure <- rbind.fill(tbl.omf, data.Procedure)
			data.Procedure <- rbind.fill(tbl.nte, data.Procedure)
			data.Procedure <- rbind.fill(tbl.f1, data.Procedure)
		#if(exists("tbl.f2"))  data.Procedure <- rbind.fill(data.Procedure, tbl.f2)
		if(exists("tbl.nte_f1")) data.Procedure <- rbind.fill(data.Procedure, tbl.nte_f1)  
		
		data.Procedure <- merge(data.Procedure, tbl.pt[, c("PatientID", "dxyear")]) 
	      
	 	df <- data.Procedure
	  	unique.dxyear<- unique(df$dxyear)
	  	unique.side<- unique(df$side)
	  	unique.site <- unique(df$site)
	  	unique.surgery_name <- unique(df$surgery_name)	  	
		unique.date<- unique(df$date)
		unique.date_locoregional<- unique(df$date_locoregional)  
		unique.date_metastatic<- unique(df$date_metastatic)

		result = list(unique.dxyear=unique.dxyear, 
	                unique.side=unique.side,
	                unique.site=unique.site,
	                unique.surgery_name=unique.surgery_name,
	                unique.date=unique.date,
					unique.date_locoregional=unique.date_locoregional,
					unique.date_metastatic=unique.date_metastatic)
	 	print(study_name)
		return(result)
	}
  #--------------------------------------------------------------------------------------------------------------------
	    df$date <- format(as.Date(df$dxyear,"%m/%d/%Y") + as.integer(df$date), "%m/%d/%Y")
	    df$date_locoregional <- format(as.Date(df$dxyear,"%m/%d/%Y") + as.integer(df$date_locoregional), "%m/%d/%Y")
	    df$date_metastatic <- format(as.Date(df$dxyear,"%m/%d/%Y") + as.integer(df$date_metastatic), "%m/%d/%Y")
#----------------------   Pathology functions Start Here   ----------------------
  Pathology.unique.request <- function(study_name){
		tbl.pt 
			 'bcr_patient_barcode' = list(name = "PatientID", data = "tcgaId"),
			 'initial_pathologic_dx_year' = list(name = "dxyear", data = "tcgaDate"), 
			 'days_to_initial_pathologic_diagnosis'  = list(name = "date", data = "upperCharacter"), #date
			 'tumor_tissue_site' = list(name = "pathDisease", data = "upperCharacter"),  
			 'histological_type'= list(name = "pathHistology", data = "upperCharacter"), 
			 'prospective_collection'= list(name = "prospective_collection", data = "upperCharacter"),
			 'retrospective_collection'= list(name = "retrospective_collection", data = "upperCharacter"), 
			 'method_initial_path_dx' = list(name = "pathMethod", data = "upperCharacter"),
			 'ajcc_tumor_pathologic_pt' = list(name = "T.Stage", data = "upperCharacter"),
			 'ajcc_nodes_pathologic_pn' = list(name = "N.Stage", data = "upperCharacter"),
			 'ajcc_metastasis_pathologic_pm' = list(name = "M.Stage", data = "upperCharacter"),
			 'ajcc_pathologic_tumor_stage'= list(name = "S.Stage", data = "upperCharacter"),
			 'ajcc_staging_edition' = list(name = "staging.System", data = "upperCharacter"),
			 'tumor_grade' = list(name = "grade", data = "upperCharacter")
		tbl.omf
			 'bcr_patient_barcode' = list(name = "PatientID", data = "tcgaId"),
			 'other_malignancy_anatomic_site' = list(name = "pathDisease", data = "upperCharacter"), 
			 'days_to_other_malignancy_dx' = list(name = "date_other_malignancy", data = "upperCharacter"), #date
			 'other_malignancy_histological_type' = list(name = "pathHistology", data = "upperCharacter"),
			 'other_malignancy_histological_type_text' = list(name = "pathHistology", data = "upperCharacter")

		data.Pathology <- rbind.fill(tbl.pt[,-match("dxyear", colnames(tbl.pt))], tbl.omf)
		data.Pathology <- merge(data.Pathology, tbl.pt[,c("PatientID", "dxyear")])
		if(any(duplicated(data.Pathology))){
		  data.Pathology <- data.Pathology[-which(duplicated(data.Pathology)), ]
		}

		 result = list(unique.dxyear=unique.dxyear,
		               unique.pathDisease=unique.pathDisease, 
                 	   unique.pathHistology=unique.pathHistology,
                 	   unique.prospective_collection=unique.prospective_collection,
                 	   unique.retrospective_collection=unique.retrospective_collection,
                 	   unique.pathMethod=unique.pathMethod,
                 	   unique.T.Stage=unique.T.Stage,
                       unique.N.Stage=unique.N.Stage,
                 	   unique.M.Stage=unique.M.Stage,
                       unique.S.Stage=unique.S.Stage,
                       unique.staging.System=unique.staging.System,
                       unique.grade=unique.grade,
                       unique.date=unique.date,
                       unique.date_other_malignancy=unique.date_other_malignancy)         
               			
               print(study_name)
  			return(result)
  }
  #-------------------------------------------------------------------------------------------------------------------------
    	df$date <- format(as.Date(df$dxyear,"%m/%d/%Y") + as.integer(df$date), "%m/%d/%Y")
    	df$date_other_malignancy <- format(as.Date(df$dxyear,"%m/%d/%Y") + as.integer(df$date_other_malignancy), "%m/%d/%Y")
#----------------------   Absent functions Start Here   -------------------------
	Absent.unique.request <- function(study_name){
	  	tbl.pt 
			 'bcr_patient_barcode' = list(name = "PatientID", data = "tcgaId"),
			 'initial_pathologic_dx_year' = list(name = "dxyear", data = "tcgaDate"),
			 'pulmonary_function_test_indicator' = list(name = "pulInd", data = "upperCharacter")
		tbl.omf 
			 'bcr_patient_barcode' = list(name = "PatientID", data = "tcgaId"),
			 'days_to_other_malignancy_dx' = list(name = "omfdx", data = "upperCharacter"),
			 'radiation_tx_indicator' = list(name = "radInd", data = "upperCharacter"),
			 'drug_tx_indicator' = list(name = "drugInd", data = "upperCharacter")
		tbl.nte 
			 'bcr_patient_barcode' = list(name = "PatientID", data = "tcgaId"),
			 'days_to_new_tumor_event_after_initial_treatment' = list(name = "omfdx", data = "upperCharacter"),
			 'new_tumor_event_dx_days_to'  = list(name = "omfdx", data = "upperCharacter"),
			 'additional_radiation_therapy' = list(name = "radInd", data = "upperCharacter"),
			 'new_tumor_event_radiation_tx' = list(name = "radInd", data = "upperCharacter"),
			 'additional_pharmaceutical_therapy' = list(name = "drugInd", data = "upperCharacter"),
			 'new_tumor_event_pharmaceutical_tx' = list(name = "drugInd", data = "upperCharacter")
		tbl.f1 
			 'bcr_patient_barcode' = list(name = "PatientID", data = "tcgaId"),
			 'new_tumor_event_dx_days_to' = list(name = "omfdx", data = "upperCharacter"),
			 'new_tumor_event_radiation_tx' = list(name = "radInd", data = "upperCharacter"),
			 'new_tumor_event_pharmaceutical_tx' = list(name = "drugInd", data = "upperCharacter")
		tbl.f2 
			 'bcr_patient_barcode' = list(name = "PatientID", data = "tcgaId"),
			 'new_tumor_event_dx_days_to' = list(name = "omfdx", data = "upperCharacter"),
			 'new_tumor_event_radiation_tx' = list(name = "radInd", data = "upperCharacter"),
			 'new_tumor_event_pharmaceutical_tx' = list(name = "drugInd", data = "upperCharacter")
		tbl.f3 
			 'bcr_patient_barcode' = list(name = "PatientID", data = "tcgaId"),
			 'new_tumor_event_dx_days_to' = list(name = "omfdx", data = "upperCharacter"),
			 'new_tumor_event_radiation_tx' = list(name = "radInd", data = "upperCharacter"),
			 'new_tumor_event_pharmaceutical_tx' = list(name = "drugInd", data = "upperCharacter")
		tbl.nte_f1 
			 'bcr_patient_barcode' = list(name = "PatientID", data = "tcgaId"),
			 'days_to_new_tumor_event_after_initial_treatment' = list(name = "omfdx", data = "upperCharacter"),
			 'new_tumor_event_dx_days_to'  = list(name = "omfdx", data = "upperCharacter"),
			 'additional_radiation_therapy' = list(name = "radInd", data = "upperCharacter"),
			 'new_tumor_event_radiation_tx' = list(name = "radInd", data = "upperCharacter"),
			 'additional_pharmaceutical_therapy' = list(name = "drugInd", data = "upperCharacter"),
			 'new_tumor_event_pharmaceutical_tx' = list(name = "drugInd", data = "upperCharacter")

	    tbl.f <- data.frame()
	    if("pulInd" %in% colnames(tbl.f)) tbl.f <- rbind.fill(tbl.f, tbl.f[,c("PatientID", "pulInd")])
	    if(exists("tbl.omf")) tbl.f <- rbind.fill(tbl.f, tbl.omf)   
	    if(exists("tbl.nte")) tbl.f <- rbind.fill(tbl.f, tbl.nte)
	    if(exists("tbl.f1")) tbl.f <- rbind.fill(tbl.f, tbl.f1)
	    if(exists("tbl.f2")) tbl.f <- rbind.fill(tbl.f, tbl.f2)
	    if(exists("tbl.f3")) tbl.f <- rbind.fill(tbl.f, tbl.f3)
	    if(exists("tbl.nte_f1")) tbl.f <- rbind.fill(tbl.f, tbl.nte_f1)
	    if(nrow(tbl.f) == 0){
	    		result = list(unique.omfdx=NA, unique.radInd=NA,
		   		   		  unique.drugInd=NA, unique.pulInd=NA)
	    }else{
	    	df <- merge(tbl.pt[,c("PatientID", "dxyear"),], tbl.f)
			unique.omfdx <- unique(df$omfdx)
			unique.radInd <- unique(df$radInd)
			unique.drugInd <- unique(df$drugInd)
			if("pulInd" %in% colnames(df)){
				unique.pulInd <- unique(df$pulInd)
			}else{
				unique.pulInd <- NA
			}
			
		   	result = list(unique.omfdx=unique.omfdx, unique.radInd=unique.radInd,
		   		   		  unique.drugInd=unique.drugInd, unique.pulInd=unique.pulInd)
	    }   
	    return(result)
	}
	#--------------------------------------------------------------------------------
		df$date <- format(as.Date(df$dxyear,"%m/%d/%Y") + as.integer(df$omfdx), "%m/%d/%Y")
#----------------------   Tests functions Start Here   --------------------------
	Tests.unique.request <- function(study_name){
	  	tbl.pt
		tbl.nte_f1 
		tbl.f1 
		tbl.f2
		tbl.f3 
		tbl.nte 
			 'bcr_patient_barcode' = list(name = "PatientID", data = "tcgaId"),
			 'initial_pathologic_dx_year' = list(name = "dxyear", data = "tcgaDate"),
			 'days_to_psa_most_recent' = list(name = "psaDate", data = "upperCharacter"),
			 'days_to_bone_scan' = list(name = "boneScanDate", data = "upperCharacter"),
			 'days_to_ct_scan_ab_pelvis' = list(name = "ctAbPelDate", data = "upperCharacter"),
			 'days_to_mri' = list(name = "mriDate", data = "upperCharacter"),
			 'idh1_mutation_test_method' =  list(name = "idh1Method", data = "upperCharacter"),
			 'idh1_mutation_found' = list(name = "idh1Found", data = "upperCharacter"),
			 'IHC' = list(name = "ihc", data = "upperCharacter"),
			 'kras_mutation_found' = list(name = "krasInd", data = "upperCharacter"),
			 'kras_mutation_identified_type' = list(name = "krasType", data = "upperCharacter"),
			 'egfr_mutation_status' = list(name = "egfrStatus", data = "upperCharacter"),
			 'egfr_mutation_identified_type' = list(name = "egfrType", data = "upperCharacter"),
			 'egfr_amplification_status' = list(name = "egfrAmp", data = "upperCharacter"),
			 'pulmonary_function_test_indicator' = list(name = "pulInd", data = "upperCharacter"),
			 'eml4_alk_translocation_status' = list(name = "elm4AlkStatus", data = "upperCharacter"),
			 'eml4_alk_translocation_variant' = list(name = "elm4AlkVar", data = "upperCharacter"),
			 'kras_mutation_codon' = list(name = "krasCodon", data = "upperCharacter"),
			 'braf_gene_analysis_indicator' = list(name = "brafInd", data = "upperCharacter"),
			 'braf_gene_analysis_result' = list(name = "brafRes", data = "upperCharacter"),
			 'cea_level_pretreatment' = list(name = "ceaTx", data = "upperCharacter"),
			 'loci_tested_count' = list(name = "lociTestCount", data = "upperCharacter"),
			 'loci_abnormal_count' = list(name = "lociAbnormalCount", data = "upperCharacter"),
			 'mismatch_rep_proteins_tested_by_ihc' = list(name = "mismatchProteinTestIhc", data = "upperCharacter"),
			 'mismatch_rep_proteins_loss_ihc' = list(name = "mismatchProteinLossIhc", data = "upperCharacter"),
			 'hpv_status_p16' = list(name = "hpvP16", data = "upperCharacter"),
			 'hpv_status_ish' = list(name = "hpvIsh", data = "upperCharacter"),
			 'psa_most_recent_results' = list(name = "psaRes", data = "upperCharacter"),
			 'bone_scan_results' = list(name = "boneScaneRes", data = "upperCharacter"),
			 'ct_scan_ab_pelvis_results' = list(name = "ctAbPelRes", data = "upperCharacter"),
			 'mri_results' = list(name = "mriRes", data = "upperCharacter"),
			 'her2_copy_number' = list(name = "her2CNV", data = "upperCharacter"),
			 'her2_fish_method' = list(name = "her2FishMethod", data = "upperCharacter"),
			 'her2_fish_status' = list(name = "her2FishStatus", data = "upperCharacter"),
			 'her2_ihc_percent_positive' = list(name = "her2IhcPercentagePos", data = "upperCharacter"),
			 'her2_ihc_score' = list(name = "her2IhcScore", data = "upperCharacter"),
			 'her2_positivity_method_text' = list(name = "her2PosMethod", data = "upperCharacter"),
			 'her2_positivity_scale_other' = list(name = "her2PosScaleOther", data = "upperCharacter"),
			 'her2_status_by_ihc' = list(name = "her2StatusIhc", data = "upperCharacter"),
			 'nte_her2_fish_define_method' = list(name = "nteHer2FishMethod", data = "upperCharacter"),
			 'nte_her2_fish_status' = list(name = "nteHer2FishStatus", data = "upperCharacter"),
			 'nte_her2_positivity_ihc_score' = list(name = "nteHer2PosIhcScore", data = "upperCharacter"),
			 'nte_her2_positivity_method' = list(name = "nteHer2PosMethod", data = "upperCharacter"),
			 'nte_her2_positivity_other_scale' = list(name = "nteHer2PosOtherScale", data = "upperCharacter"),
			 'nte_her2_signal_number' = list(name = "nteHer2SignalNum", data = "upperCharacter"),
			 'nte_her2_status' = list(name = "nteHer2Status", data = "upperCharacter"),
			 'nte_her2_status_ihc_positive' = list(name = "nteHer2StatusIhcPos", data = "upperCharacter"),
			 'nte_er_ihc_intensity_score' = list(name = "nteEstroIhcScore", data = "upperCharacter"),
			 'nte_er_positivity_define_method' = list(name = "nteEstroPosMethod", data = "upperCharacter"),
			 'nte_er_positivity_other_scale' = list(name = "nteEstroPosOtherScale", data = "upperCharacter"),
			 'nte_er_status' = list(name = "nteEstroStatus", data = "upperCharacter"),
			 'nte_er_status_ihc_positive' = list(name = "nteEstroStatusIhcPos", data = "upperCharacter"),
			 'nte_pr_ihc_intensity_score' = list(name = "nteProgIhcScore", data = "upperCharacter"),
			 'nte_pr_positivity_define_method' = list(name = "nteProgPosMethod", data = "upperCharacter"),
			 'nte_pr_positivity_other_scale' = list(name = "nteProgPosOtherScale", data = "upperCharacter"),
			 'nte_pr_status_by_ihc' = list(name = "nteProgStatusIhc", data = "upperCharacter"),
			 'nte_pr_status_ihc_positive' = list(name = "nteProgStatusIhcPos", data = "upperCharacter"),
			 'pr_positivity_define_method' = list(name = "ProgPosMethod", data = "upperCharacter"), 
			 'pr_positivity_ihc_intensity_score' = list(name = "ProgPosIhcScore", data = "upperCharacter"),
			 'pr_positivity_scale_other' = list(name = "ProgPosScaleOther", data = "upperCharacter"),
			 'pr_positivity_scale_used' = list(name = "ProgPosScaleUsed", data = "upperCharacter"),
			 'pr_status_by_ihc' = list(name = "ProgStatusIhc", data = "upperCharacter"),
			 'pr_status_ihc_percent_positiv' = list(name = "ProgStatusIhcPercentagePos", data = "upperCharacter"),
			 'her2_and_cent17_cells_count' = list(name = "her2Cent17CellsCount", data = "upperCharacter"),
			 'her2_and_cent17_scale_other' = list(name = "her2Cent17ScaleOther", data = "upperCharacter"),
			 'her2_cent17_counted_cells_count' = list(name = "her2Cent17CountedCellsCount", data = "upperCharacter"),
			 'her2_cent17_ratio' = list(name = "her2Cent17Ratio", data = "upperCharacter"),
			 'nte_cent_17_her2_ratio' = list(name = "nteCent17Her2Ratio", data = "upperCharacter"),
			 'nte_cent_17_signal_number' = list(name = "nteCent17SignalNum", data = "upperCharacter"),
			 'nte_cent17_her2_other_scale' = list(name = "nteCent17Her2OtherScale", data = "upperCharacter")

	  	tbl <- list()
	  	if(length(tbl.pt) > 2){
	  		tbl <- tbl.pt[,-match("dxyear", colnames(tbl.pt))]
	  	}
	  	if(exists("tbl.f1") && length(tbl.f1) > 1) tbl <- rbind.fill(tbl, tbl.f1)
	  	if(exists("tbl.f2") && length(tbl.f2) > 1) tbl <- rbind.fill(tbl, tbl.f2)
	  	if(exists("tbl.f3") && length(tbl.f3) > 1) tbl <- rbind.fill(tbl, tbl.f3)
	  	if(exists("tbl.nte") && length(tbl.nte) > 1) tbl <- rbind.fill(tbl, tbl.nte)
	  	if(exists("tbl.nte_f1") && length(tbl.nte_f1) > 1) tbl <- rbind.fill(tbl, tbl.nte_f1)
	    if(length(tbl) == 0){
	    	print(c(study_name, length(tbl)))
	    	return(tbl)
	    }else{	
			unique.Result <- unique(unlist(tbl[,-1]))
		   	result = list(unique.Result=unique.Result)
		  	return(result)
		}
	}
	#--------------------------------------------------------------------------------
	Tests.mapping.type <- function(df, df.methods){
		if(exists("type")){
			rm(type)		
		}
		type <- NA
		if(missing(df.methods)){
			if(length(grep("psa", colnames(df[2]), ignore.case=TRUE)) > 0 ) {type <- "PSA"}
			if(length(grep("boneScan", colnames(df[2]), ignore.case=TRUE)) > 0 ){type <- "BONE SCAN"}
			if(length(grep("ctAbPel", colnames(df[2]), ignore.case=TRUE)) > 0 ){type <- "CT SCAN"}
			if(length(grep("mri", colnames(df[2]), ignore.case=TRUE)) > 0 ){type <- "MRI"}
			if(length(grep("ihc", colnames(df[2]), ignore.case=TRUE)) > 0 ){type <- "IHC"}
			if(length(grep("pul", colnames(df[2]), ignore.case=TRUE)) > 0 ){type <- "PULMONARY"}
			if(length(grep("p16", colnames(df[2]), ignore.case=TRUE)) > 0 ){type <- "P16"}
			if(length(grep("ish", colnames(df[2]), ignore.case=TRUE)) > 0 ){type <- "ISH"}
			if(length(grep("fish", colnames(df[2]), ignore.case=TRUE)) > 0 ){type <- "FISH"}
			if(length(grep("cellsCount", colnames(df[2]), ignore.case=TRUE)) > 0 ){type <- "CELLS COUNT"}
			if(!exists("type")) type <- NA
			df$Type <- rep(type, nrow(df))
			return(df)
		}else{
				methods <- gsub("Method", "", colnames(df.methods))
				for(i in 1:length(methods)){
					if(length(grep(methods[i], colnames(df)[2], ignore.case=TRUE)) > 0){
						type <- df.methods[,i]
					}
				}
				if(is.na(type)){
					if(length(grep("psa", colnames(df[2]), ignore.case=TRUE)) > 0 ) {type <- rep("PSA", nrow(df))}
					if(length(grep("boneScan", colnames(df[2]), ignore.case=TRUE)) > 0 ){type <- rep("BONE SCAN", nrow(df))}
					if(length(grep("ctAbPel", colnames(df[2]), ignore.case=TRUE)) > 0 ){type <- rep("CT SCAN", nrow(df))}
					if(length(grep("mri", colnames(df[2]), ignore.case=TRUE)) > 0 ){type <- rep("MRI", nrow(df))}
					if(length(grep("ihc", colnames(df[2]), ignore.case=TRUE)) > 0 ){type <- rep("IHC", nrow(df))}
					if(length(grep("pul", colnames(df[2]), ignore.case=TRUE)) > 0 ){type <- rep("PULMONARY", nrow(df))}
					if(length(grep("p16", colnames(df[2]), ignore.case=TRUE)) > 0 ){type <- rep("P16", nrow(df))}
					if(length(grep("ish", colnames(df[2]), ignore.case=TRUE)) > 0 ){type <- rep("ISH", nrow(df))}
					if(length(grep("fish", colnames(df[2]), ignore.case=TRUE)) > 0 ){type <- rep("FISH", nrow(df))}
					if(length(grep("cellsCount", colnames(df[2]), ignore.case=TRUE)) > 0 ){type <- rep("CELLS COUNT", nrow(df))}
				}
				df$Type <- type
				return(df)
		}
	}
	#--------------------------------------------------------------------------------
	Tests.mapping.testAndDate <- function(df, df.dates){
		testNames <- c("KRAS", "EGFR", "Pulmonary_function", "IDH1", "EML4_ALK", "BRAF", "CEA", "LOCI", "MISMATCHED_PROTEIN",
  					"HPV_P16", "HPV_ISH", "PSA", "BONE_SCAN", "CT_SCAN_AB_PELVIS", "MRI", "HER2", "ESTROGEN", 
  					"PROGESTERONE_RECEPTOR", "CENTROMERE_17")
	    searchKeyWords <- c("kras", "egfr", "pul", "idh1", "elm4Alk", "braf", "cea", "loci", "mismatchProtein", 
	    					"hpvP16", "hpvIsh", "psa", "boneScan", "ctAbPel", "mri", "her2", "estro", "prog", "cent17")
	    test <- NA
	    key <- NA
	    for(i in 1:length(searchKeyWords)) {
	    	if(length(grep(searchKeyWords[i], colnames(df[2]),ignore.case=TRUE)) > 0){
	    		key <- searchKeyWords[i]
	    		test <- testNames[i]
	    	}
	    }
	    if(!missing(df.dates)) {
			if(length(grep(key, colnames(df.dates))) > 0) {
	    		df$Date <- df.dates[,grep(key, colnames(df.dates))]
	    		df$Test <- rep(test, nrow(df))
				return(df)
			}
		}else{
			df$Test <- rep(test, nrow(df))
			return(df)
		}
	}
	#--------------------------------------------------------------------------------
	Tests.mapping.date.Calculation <- function(df){
		df$Date <- format(as.Date(df$dxyear,"%m/%d/%Y") + as.integer(df$Date), "%m/%d/%Y")
		return(df)
	}
	#--------------------------------------------------------------------------------	
} # End of Absent Native Functions
################################################     Step 4: Generate Result    ###########################################
create.DOB.records <- function(study_name, ptID){
    				return(list(PatientID=PatientID, PtNum=PtNum, study=study_name, Name="Birth", 
    				 			Fields=list(date=date, gender=gender, race=race, ethnicity=ethnicity)))	
    				})
#--------------------------------------------------------------------------------------------------------------------------
create.Diagnosis.records <- function(study_name, ptID){
    				return(list(PatientID=PatientID, PtNum=PtNum, study=study_name, Name="Diagnosis", 
    				 			Fields=list(date=date, disease=disease, siteCode=siteCode)))
    				})
#--------------------------------------------------------------------------------------------------------------------------
create.Chemo.records <- function(study_name,  ptID){
					return(list(PatientID=PatientID, PtNum=PtNum, study=study_name, Name="Drug", 
		    				 			Fields=list(date=date, agent=agent, therapyType=therapyType, intent=intent,
		    				 				        dose=dose, units=units, totalDose=totalDose, totalDoseUnits=totalDoseUnits,
		    				 				        route=route,cycle=cycle)))
		    				})
#--------------------------------------------------------------------------------------------------------------------------
create.Rad.records <- function(study_name,  ptID){
		    				return(list(PatientID=PatientID, PtNum=PtNum, study=study_name, Name="Radiation", 
		    				 			Fields=list(date=date, therapyType=radType, intent=intent, 
		    				 						target=target, totalDose=totalDose, totalDoseUnits=totalDoseUnits, 
		    				 						numFractions=numFractions)))
		    				})
#--------------------------------------------------------------------------------------------------------------------------
create.Status.records <- function(study_name,  ptID){
 	for(i in 1:nrow(tbl.pt)){
 		tmpDF <- subset(data.Status, PatientID == tbl.pt$PatientID[i])
 		tmpDF <- tmpDF[order(as.integer(tmpDF$date), decreasing=TRUE, na.last=TRUE),]
 		if(nrow(tmpDF[which(tmpDF$date == tmpDF[1,]$date), ]) > 1){
 			tmpDup <- tmpDF[which(tmpDF$date == tmpDF[1,]$date), ]
 			tmpDF[1, "vital"] 		= 	ifelse(any(duplicated(tmpDup[,"vital"])), tmpDup[1, "vital"], paste(tmpDup[, "vital"]))
 			tmpDF[1, "tumorStatus"] = 	ifelse(any(duplicated(tmpDup[,"tumorStatus"])), tmpDup[1, "tumorStatus"], paste(tmpDup[, "tumorStatus"], collapse=";"))
		}
		recentTbl <- rbind.fill(recentTbl, tmpDF[1,])
 	}

 	data.Status <- Status.mapping.date.Calculation(recentTbl)
 	data.Status[order(data.Status$PatientID, data.Status$date, data.Status$vital, data.Status$tumorStatus),] -> data.Status

    				return(list(PatientID=PatientID, PtNum=PtNum, study=study_name, Name="Status", 
    				 			Fields=list(date=date, status=vital, tumorStatus=tumorStatus)))
    				})
#--------------------------------------------------------------------------------------------------------------------------
create.Progression.records <- function(study_name,  ptID){
	    				return(list(PatientID=PatientID, PtNum=PtNum, study=study_name, Name="Progression", 
	    				 			Fields=list(date=date, event=event, number=number)))
	    				})
#--------------------------------------------------------------------------------------------------------------------------
create.Absent.records <- function(study_name,  ptID){
	    				return(list(PatientID=PatientID, PtNum=PtNum, study=study_name, Name="Absent", 
	    				 			Fields=list(date=date, Radiation=rad, Drug=drug, Pulmonary=pul)))
	    				})
#--------------------------------------------------------------------------------------------------------------------------
create.Tests.records <- function(study_name,  ptID){
	    				return(list(PatientID=PatientID, PtNum=PtNum, study=study_name, Name="Tests", 
	    				 			Fields=list(date=date, Type=type, Test=test, Result=result)))
	    				})
}
#--------------------------------------------------------------------------------------------------------------------------
create.Encounter.records <- function(study_name,  ptID){
			    return(list(PatientID=PatientID, PtNum=PtNum, study=study_name, Name="Encounter", 
			                Fields=list(encType=encType, KPS=KPS, ECOG=ECOG, height=height,
			                            weight=weight, prefev1.ratio=prefev1.ratio, prefev1.percent=prefev1.percent, postfev1.ratio=postfev1.ratio,
			                            postfev1.percent=postfev1.percent,carbon.monoxide.diffusion=carbon.monoxide.diffusion)))
			    })
#--------------------------------------------------------------------------------------------------------------------------
create.Procedure.records <- function(study_name,  ptID){
		      			return(list(PatientID=PatientID, PtNum=PtNum, study=study_name, Name="Procedure", 
		                  			Fields=list(date=date,date_locoregional=date_locoregional,date_metastatic=date_metastatic,name=name,site=site,side=side)))
#--------------------------------------------------------------------------------------------------------------------------  
create.Pathology.records <- function(study_name,  ptID){
		    return(list(PatientID=PatientID, PtNum=PtNum, study=study_name, Name="Pathology", 
		                Fields=list(pathDisease=pathDisease,
		                            pathHistology=pathHistology,
		                            prospective_collection=prospective_collection,
		                            retrospective_collection = retrospective_collection,
		                            pathMethod=pathMethod,
		                            T.Stage=T.Stage,
		                            N.Stage=N.Stage,
		                            M.Stage=M.Stage,
		                            S.Stage=S.Stage,
		                            staging.System=staging.System,
		                            grade=grade,
		                            date=date,
									date_other_malignancy=date_other_malignancy)))            
		  })
######################################    Step 5: Generate Result By Organ Site   #########################################
create.STUDY.records <- function(study_name){
	dob.events <- create.DOB.records(study_name)
	diagnosis.events <- create.Diagnosis.records(study_name)
	chemo.events <- create.Chemo.records(study_name)
	radiation.events <- create.Rad.records(study_name)
	status.events <- create.Status.records(study_name)
	progression.events <- create.Progression.records(study_name)
	absent.events <- create.Absent.records(study_name)
	tests.events <- create.Tests.records(study_name)
	encounter.events <- create.Encounter.records(study_name)
	procedure.events <- create.Procedure.records(study_name)
	pathology.events <- create.Pathology.records(study_name)

	events <- append(dob.events, chemo.events)
    events <- append(events, diagnosis.events)
    events <- append(events, status.events)
    events <- append(events, progression.events)
    events <- append(events, radiation.events)
    events <- append(events, procedure.events)
    events <- append(events, encounter.events)
    events <- append(events, pathology.events)
    events <- append(events, absent.events)
    events <- append(events, tests.events)
    print(table(unlist(lapply(events, function(e) e["Name"]))))
    events
}

brca <- create.STUDY.records(studies[1])
coad <- create.STUDY.records(studies[2])
gbm  <- create.STUDY.records(studies[3])
hnsc <- create.STUDY.records(studies[4])
lgg  <- create.STUDY.records(studies[5])
luad <- create.STUDY.records(studies[6])
lusc <- create.STUDY.records(studies[7])
prad <- create.STUDY.records(studies[8])
read <- create.STUDY.records(studies[9])
sarc <- create.STUDY.records(studies[10])
laml <- create.STUDY.records(studies[11])
blca <- create.STUDY.records(studies[12])
paad <- create.STUDY.records(studies[13])

# run through all studies by Feature
lapply(studies, create.DOB.records)
lapply(studies, create.Diagnosis.records)
lapply(studies, create.Chemo.records) #equipped with filling func
lapply(studies, create.Rad.records) #equipped with filling func
lapply(studies, create.Status.records) #equipped with filling func
lapply(studies, create.Progression.records)
lapply(studies, create.Absent.records)
lapply(studies, create.Procedure.records)
lapply(studies, create.Encounter.records)
lapply(studies, create.Pathology.records) 
lapply(studies, create.Tests.records) 

###########################################    Step 6: UnitTests By Feature  ###############################################
source("TCGAPipelingRefactoring_testingSuite.R")