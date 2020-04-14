const { Drive } = require('../models/drive');

async function postDriveFile(req,res){
    // const drive =  await Drive.find()
    const drive = new Drive({
        name: 'initial',
        reference_primary: 'ref_one.txt',
        reference_sec:'ref_two.txt',
        createdAt: '15555543',
        graphs: false,
        device_primary: 'IMU 300RI',
        device_secondary:'Novatel'
    })
    
    drive.save(function(err){
        if (err) return err;
        res.status(200).json({
            success:true
        })
    })
}

module.exports = {
    postDriveFile
}