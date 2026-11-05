cd /tmp/claude-0/wt239
source /etc/profile.d/planetgen-test.sh
export PYTHONPATH=$PWD/src PLANETGEN_BAND_TABLE_CACHE=/tmp/claude-0/perf233-bandcache PLANETGEN_TILE_CACHE_DIR=/tmp/claude-0/perf233/tiles PLANETGEN_CACHE_MAX_MB=100000 PLANETGEN_MYSQL_DATABASE_PREFIX=perf
date +%T
python /tmp/claude-0/perf233/gen.py perf233 warm --progress 2>&1 | grep -E "^Phase|^Cache|^error|Traceback|Error"
date +%T
