CLASSPATH=.
for i in *.jar; do
        CLASSPATH=$CLASSPATH:$i
done

CLASSPATH=$CLASSPATH:bin
export CLASSPATH

java -Xss1g -Xmx10g -XX:-UseGCOverheadLimit mining.EPM $1 $2 $3 $4

