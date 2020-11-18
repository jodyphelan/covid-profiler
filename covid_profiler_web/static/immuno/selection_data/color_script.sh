cat selection_data.raw | while read -r line

do

echo $line | awk -F "," '{

if ($8 >=0.1 && $8 <0.16)
  print $0",#95676a"

  if ($8 >=0.16 && $8 <0.32)
    print $0",#aa5c55"

    if ($8 >=0.32 && $8 <0.48)
      print $0",#bf5040"

      if ($8 >=0.48 && $8 <0.64)
        print $0",#d4442b"

        if ($8 >=0.64 && $8 <0.8)
          print $0",#ea3915"

          if ($8 >=0.8 && $8 <1)
          	print $0",#ff2d00"

            if ($8 >=0.00000000000001 && $8 <0.1)
              print $0",#807380"

              if ($8 ==0)
                print $0",#807380"

              if ($8 <=-0.00000000000001 && $8 >-0.1)
                print $0",#807380"

            if ($8 <= -0.1 && $8 > -0.16)
              print $0",#6a7f95"

              if ($8 <= -0.16 && $8 > -0.32)
                print $0",#558aaa"

                if ($8 <= -0.32 && $8 > -0.48)
                  print $0",#4096bf"

                  if ($8 <= -0.48 && $8 > -0.64)
                    print $0",#2ba2d4"

                    if ($8 <= -0.64 && $8 > -0.8)
                      print $0",#15adea"

                      if ($8 <= -0.8 && $8 > -1)
                      	print $0",#00b9ff"

}'



done > parsed_color_mutations.tmp

cat ../gallery05_SNP01.js | awk -F '"' '{print $2$4$0}' |  grep chr | sort > SNP_delim.tmp

cat parsed_color_mutations.tmp | sed 's/,,/,/g'  |  awk -F "," '{print $1$2"\t"$8"\""$9}' | sort > SEL_delim.tmp

join -1 1  -2 1 -a1  SNP_delim.tmp SEL_delim.tmp
