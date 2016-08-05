const connectionString = 'mongodb://localhost:27017/oncoscape';
const dataDir = "../data";
const comongo = require('co-mongodb')
var co = require('co');
const _ = require('underscore');
const assert = require('assert');
const fs = require('fs');
const rl = require('readline');
const request = require('sync-request');
const xml = require('xmldom');
const execSync = require('child_process').execSync;
var collection, data, i, j;


var onerror = function(e){
	console.log("SHIT!");
	console.log(e);
}


co(function *() {
	
	// db is just a regular Db instance from the native driver. 
	db = yield comongo.client.connect(connectionString);

	var d = JSON.parse(fs.readFileSync("oncoscape.brca_mut_cBio_wxs"));
	collection = yield comongo.db.collection(db, "oncoscape.brca_mut_cBio_wxs");
	for (var i=0; i<d.rows.length; i++){
		console.log(d.rows[i])
		var o = {gene:(d.rows[i]!=null) ? d.rows[i].split("|")[0] : null, data:d.data[i]};
		yield collection.insert(o, {w:"majority"});	
	}

	
	yield collection.insert(d1, {w:"majority"});
	
	 
	yield comongo.db.close(db);
}).catch(onerror);