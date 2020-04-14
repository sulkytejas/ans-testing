const { Drive } = require('../models/drive');

async function getAllFiles(req,res){
    const drive =  await Drive.find()
    let files = drive.map(f => (
        {
            name:f.name,
            // createdAt:f.createdAt,
            reference_primary:f.reference_primary,
            reference_sec:f.reference_sec,
            devicePrimary: f.device_primary,
            deviceSecondary: f.device_secondary,
            graphs: f.graphs
        }
    ))
    
    res.json({
        success: true,
        files
    })
}

module.exports = {
    getAllFiles
}