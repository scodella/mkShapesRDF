#!/bin/sh -x
git checkout master
git fetch
#git pull
git checkout run3base
git merge origin/master

