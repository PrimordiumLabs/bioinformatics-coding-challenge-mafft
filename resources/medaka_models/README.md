This folder is to store medaka models we want to use. They are downloaded from https://github.com/nanoporetech/medaka/tree/master/medaka/data

This ensures that when the ONT server goes down our pipeline does not start failing. 

The models ending with "_pt" are for medaka version 2 (they are the pytorch versions)

If you want to switch to a different medaka model, just add it here and change the medaka_model in config.yaml. 