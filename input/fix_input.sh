for i in ./*/*.dat; do grep -v -e 'N =' -e 'Re =' -e 'We =' -e 'beta =' $i > $i.temp;mv $i.temp $i; done
