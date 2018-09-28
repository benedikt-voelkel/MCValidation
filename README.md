# Validation package for VMC

This package puts together a few things for validation of MC engines run with ROOT's Virtual Monte Carlo (VMC) package. For evaluation the [MCStepLogger](https://github.com/AliceO2Group/AliceO2/tree/dev/Utilities/MCStepLogger) developed in the ALICE O2 software framework is used. Since it is inherently independent of O2, it was extracted and adapted and the version used can be found [here](https://github.com/benedikt-voelkel/MCStepLogger).

## Requirements

The minimal requirements are

* [Python](https://www.python.org/downloads/release/python-2715/), version 2.7.X together with Python's package manager `pip`.
* `git`
Everything else is taken care of by the ALICE build tool [aliBuild](https://github.com/alisw/alibuild). Note, that everything will happen in the top directory of your cloned repository.

[ROOT](https://github.com/root-project/root), preferably a version >= 6.X (tested with 6.12.06),
* [Boost](https://www.boost.org/users/download/), version >= 1.59.

## Setup a validation
In the top directory there is a `setup.sh` which you need to run. If not present, it will clone the `git` repositories necessary for all following steps and it will locally install `aliBuild` for you. Also, it will export some environment variables

All sub-directories contain a `setup.sh` and `CMakeLists.txt` to easily setup and build the project. In the following the validation scenarios are explained in detail.

## Moving tracks between stacks
So far, there is one main validation workflow in [stackMove](./stackMove). The VMC package is further developed to allow for different engines running concurrently. In this scenario it is then possible to distribute the simulation to different engines based on arbitrary criteria. The particle transport through the entire geometry setup is then done in one run by making sure that tracks are properly stopped and transferred to another engine as soon as certain user criteria are met.

In order to investigate the performance of the code managing the stopping, criteria evaluation and transferring of the respective tracks, the following is done:
Choose a transport engine (here either GEANT3 or GEANT4):
* Run the engine for a given geometry and monitor performance.
* Compare performance to the setup where the same engine is used but where at certain points a track is stopped, transferred to the same engine and started again.

same engine is run standalone and then compared to a setup where the same engine  comparing the performance of running a single engine vs. this engine in concurrent mode
