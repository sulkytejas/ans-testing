const { Router } = require ('express');
const { getAllFiles } =  require ('../controllers/getAllFiles');
const { postDriveFile } =  require ('../controllers/postDriveFile');
const { deleteDriveFile } =  require ('../controllers/deleteDriveFile');

const driveFileRouter = Router();

driveFileRouter.get('/api/drive/logs', getAllFiles);
driveFileRouter.post('/api/drive/insert', postDriveFile);
driveFileRouter.delete('/api/drive/logs/:name/delete', deleteDriveFile);

module.exports = {
    driveFileRouter
}

