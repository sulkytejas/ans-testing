const express = require ('express');
const app = express();
const { BlobServiceClient,ContainerClient } = require('@azure/storage-blob');
const cors = require('cors');
const mongoose = require('mongoose');
const {spawn} = require('child_process');

const AZURE_STORAGE_CONNECTION_STRING = process.env.AZURE_STORAGE_CONNECTION_STRING;

const { Drive } = require('./models/drive');

require('dotenv').config();
app.use(express.json());
app.use(express.urlencoded({
  extended: true
}));
app.use(cors());

mongoose.Promise = global.Promise;
mongoose.connect(process.env.DATABASE, {
  auth: {
    user: process.env.COSMODDB_USER,
    password: process.env.COSMOSDB_PASSWORD
  },
  useCreateIndex: true,
  useUnifiedTopology: true,
  useNewUrlParser: true
})
.then(() => {
  // realtimeBaseObservation();
  console.log(`Database sucessfull connected`)
})
.catch((err) => console.error(err));



/* Api call to fetch file from azzure */
const  { driveFileRouter } = require ('./routes/driveFileRouter');
let complete = false

app.use(driveFileRouter)

 /* start processing */

app.get('/api/drive/process/:name', async function (req,res){
  const file = await Drive.findOne({
    name: `${req.params.name}`
  });
  
  const python = spawn('python', ['./post_processing/txtProcessing.py',`${file.name}`]);
  
  file.graphs = true;
  await file.save()

  res.json({
    success: true
  });
});

/* trigger to finish processing */

let clients = []
app.get('/api/drive/complete/', function (req,res){
 
 clients.forEach(c => c.write(`data: completed\n\n`))
  res.json({
    success: true
  });

})



/* api triger to notify web about the file */
app.get('/added', function(req, res) {
  clients.push(res)
  res.set({
    "Content-Type": "text/event-stream",
    "Cache-Control": "no-cache",
    Connection: "keep-alive",

    // enabling CORS
    "Access-Control-Allow-Origin": "*",
    "Access-Control-Allow-Headers":
      "Origin, X-Requested-With, Content-Type, Accept",
  })

  res.write(`data:start\n\n`)


  req.on('close', function(){
    console.log('close')
    res.end();
  })

});

const port = process.env.PORT;
app.listen(port, () => console.log(`App running on ${port}`))

