mongo <- connect.to.mongo()
dataset = "brca"
dataType = "cnv"
source = "cBio"
name = "fix me!!"

lookup.ns <-  "oncoscape.lookup_oncoscape_datasources"
approach = "buffer" # list

if(approach == "list"){
  
  query <- list(disease=dataset )
  datasource <- mongo.find.all(mongo,lookup.ns, query)[[1]]

  if(length(datasource) == 0)
    datasource$disease = dataset

  ## add collection
    add.collection <- data.frame(source=source, type=dataType, collection=name)
    if("molecular" %in% names(datasource)){
      datasource$molecular	<- rbind(datasource$molecular, add.collection)
    } else { datasource$molecular <- add.collection   }
  
  ## insert into mongo
    mongo.update(mongo, lookup.ns, query, datasource)
}

if(approach=="buffer"){

  query <- list("disease"=dataset)
  datasource <- mongo.find.one(mongo, lookup.ns, query)

  if(length(datasource)==0){
    data.list <- list();
    data.list$disease = dataset
  }else{
    data.list <- mongo.bson.to.list(datasource)
  }
  add.collection <- list(data.frame(source=source, type=dataType, collection=name))
  if("molecular" %in% names(data.list)){
    data.list$molecular <- c(data.list$molecular, add.collection)
  }else{data.list$molecular <- add.collection}
  
  
  ## insert into mongo
  mongo.update(mongo, lookup.ns, query, data.list, mongo.update.upsert)
  
}