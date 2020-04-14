const { Drive } = require('../models/drive');

async function deleteDriveFile(req,res){
    const name = req.params.name
        
      const count = await Drive.count({
        name
      });
    
      if (count <= 0) {
        return res.json({
          success: false,
        });
      }

    const drive = await Drive.deleteOne({
        name
      });
     
      return res.json({
        success: true
      });
}

module.exports = {
    deleteDriveFile
}