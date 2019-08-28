# INTRO

# Basic Usage

Usage: python Intro.py -m [ Dstat | RND | Dxy ] ...

## Dstat

Usage: python Intro.py -m Dstat -v [VCF] -i [INFO] -o [OUT-POP] -t [THR-POP] -p [POP1] -q [POP2] -c [CHR] -w [WINDOW] -s [SLIDE]

![Dstat](https://user-images.githubusercontent.com/49300659/63832735-69513b80-c9ac-11e9-93fe-0b656cb363eb.png)

    Example: python Intro.py -m Dstat \
    
                         --input ExampleData/TestSet_234.chr10.vcf.gz \
                         
                         --info ExampleData/TestSet_234.sample_group_info.txt \
                         
                         --popO Wild_GroupC \
                         
                         --popT Wild_Mix_GroupB \
                         
                         --pop1 Cultivar_A \
                         
                         --pop2 Cultivar_OutGroup \
                         
                         --chr 10 \
                         
                         --window 50000 \
                         
                         --slide 25000


## RND

Usage: python Intro.py -m RND -v [VCF] -i [INFO] -o [OUT-POP] -p [POP1] -q [POP2] -c [CHR] -w [WINDOW] -s [SLIDE]

![RND](https://user-images.githubusercontent.com/49300659/63832750-71a97680-c9ac-11e9-8a63-f413eec203bf.png)

    Example: python Intro.py -m RND \

                         --input ExampleData/TestSet_234.chr10.vcf.gz \
                         
                         --info ExampleData/TestSet_234.sample_group_info.txt \
                         
                         --popO Wild_Mix_GroupB \
                         
                         --pop1 Cultivar_A \
                         
                         --pop2 Cultivar_OutGroup \
                         
                         --chr 10 \
                         
                         --window 50000 \
                         
                         --slide 25000


## Dxy

Usage: python Intro.py -m Dxy -v [VCF] -i [INFO] -p [POP1] -q [POP2] -c [CHR] -w [WINDOW] -s [SLIDE]

![Dxy](https://user-images.githubusercontent.com/49300659/63832768-7837ee00-c9ac-11e9-805f-f955aa8e5f5b.png)

    Example: python Intro.py -m Dxy \

                         --input ExampleData/TestSet_234.chr10.vcf.gz \
                         
                         --info ExampleData/TestSet_234.sample_group_info.txt \
                         
                         --pop1 Cultivar_A \
                         
                         --pop2 Cultivar_OutGroup \
                         
                         --chr 10 \
                         
                         --window 50000 \
                         
                         --slide 25000

## Plot

Usage: Rscript IntroPlot.R [METHOD] [RESULT_FILE] [CHR]

![R](https://user-images.githubusercontent.com/49300659/63834056-5ab85380-c9af-11e9-8bfe-ee1bcafbbda8.png)


# Calculation Flow
![그림1](https://user-images.githubusercontent.com/49300659/63830448-1fb22200-c9a7-11e9-86f5-ba709246719c.jpg)

# Requirement

Python 3.0 and numpy library are requiered. 


# Contact

jaeyoonkim72@gmail.com

# License

Copyright number in Korea: C-2017-024343

# Reference

Feder, J. L., Xie, X., Rull, J., Velez, S., Forbes, A., Leung, B., ... & Aluja, M. (2005). Mayr, Dobzhansky, and Bush and the complexities of sympatric speciation in Rhagoletis. Proceedings of the National Academy of Sciences, 102(suppl 1), 6573-6580.

Nei, M., & Li, W. H. (1979). Mathematical model for studying genetic variation in terms of restriction endonucleases. Proceedings of the National Academy of Sciences, 76(10), 5269-5273.

Patterson, N., Moorjani, P., Luo, Y., Mallick, S., Rohland, N., Zhan, Y., ... & Reich, D. (2012). Ancient admixture in human history. Genetics, 192(3), 1065-1093.
