#!/bin/sh -x
git checkout base
git fetch
#git pull
git checkout BTagPerf
git merge origin/base

