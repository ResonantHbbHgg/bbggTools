#!/bin/bash

watch -n 1 'ps aux | grep ${USER} | grep cmsRun | grep Make | grep -v ps'
