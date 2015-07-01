#!/bin/bash

watch -n 1 'ps aux | grep ${USER} | grep cmsRun | grep MakePlots.py | grep -v ps'
