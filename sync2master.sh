#!/bin/sh -x
git checkout master
git fetch
#git pull
git checkout base
git merge origin/master

