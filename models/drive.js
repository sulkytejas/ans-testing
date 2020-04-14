const mongoose = require('mongoose');
var ObjectId = mongoose.Schema.Types.ObjectId;

const driveSchema = mongoose.Schema({
    name:{
        type: String,
        required: true
    },
    reference_primary:{
        type: String,
        require: false
    },
    reference_sec:{
        type: String,
        require: false
    },
    createdAt:{
        type:Date,
        required: false
    },
    graphs:{
        type:Boolean,
        required: true    
    },
    device_primary:{
        type:String,
        required: false    
    },
    device_secondary:{
        type:String,
        required: false    
    }
});

const Drive = mongoose.model('drive',driveSchema,'drive');

module.exports = { Drive };