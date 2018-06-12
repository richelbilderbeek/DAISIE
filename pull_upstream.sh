#!/bin/bash
git remote add upstream git://github.com/rsetienne/DAISIE.git
git fetch origin -v; git fetch upstream -v; git merge upstream/master