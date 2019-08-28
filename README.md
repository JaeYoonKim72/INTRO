# INTRO

INTRO is a free software that can identify genomic regions introgressed between two groups.

Source code was written by Jae-Yoon Kim using Python language, and supported on windows and linux platform.

This software is based on three statistics, D-statistics, RND, and Dxy, raised by Patterson, Feder, and Nei.

The description of these three statistics is as follows.

![그림5](https://user-images.githubusercontent.com/49300659/63851385-2bb2d980-c9d1-11e9-8848-c04ab430c38b.png)


# Basic Usage

Usage: python Intro.py -m [ Dstat | RND | Dxy ] ...

## 1) Dstat

Usage: python Intro.py -m Dstat -v [VCF] -i [INFO] -o [OUT-POP] -t [THR-POP] -p [POP1] -q [POP2] -c [CHR] -w [WINDOW] -s [SLIDE]

![Dstat](https://user-images.githubusercontent.com/49300659/63832735-69513b80-c9ac-11e9-93fe-0b656cb363eb.png)

    Example: python Intro.py -m Dstat \
    
                         --input ExampleData/TestSet_234.chr10.vcf.gz \
                         
                         --info ExampleData/TestSet_234.sample_group_info.txt \
                         
                         --popO Wild_Mix_GroupA \
                         
                         --popT Cultivar_B \
                         
                         --pop1 Cultivar_A \
                         
                         --pop2 Cultivar_OutGroup \
                         
                         --chr 10 \
                         
                         --window 50000 \
                         
                         --slide 25000


## 2) RND

Usage: python Intro.py -m RND -v [VCF] -i [INFO] -o [OUT-POP] -p [POP1] -q [POP2] -c [CHR] -w [WINDOW] -s [SLIDE]

![RND](https://user-images.githubusercontent.com/49300659/63832750-71a97680-c9ac-11e9-8a63-f413eec203bf.png)

    Example: python Intro.py -m RND \

                         --input ExampleData/TestSet_234.chr10.vcf.gz \
                         
                         --info ExampleData/TestSet_234.sample_group_info.txt \
                         
                         --popO Wild_Mix_GroupA \
                         
                         --pop1 Cultivar_A \
                         
                         --pop2 Cultivar_OutGroup \
                         
                         --chr 10 \
                         
                         --window 50000 \
                         
                         --slide 25000


## 3) Dxy

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

# Plot Usage
## 1) Usage

Usage: Rscript IntroPlot.R [METHOD] [RESULT_FILE] [CHR]

![R](https://user-images.githubusercontent.com/49300659/63834107-77ed2200-c9af-11e9-889a-584d42b882fb.png)

    Example: Rscript IntroPlot.R \
                     Dstat \
                     Dtest.50000.25000.10.Cultivar_B-Cultivar_A-Cultivar_OutGroup-Wild_Mix_GroupA.txt \
                     10

    Example: Rscript IntroPlot.R \
                     RND \
                     RND.50000.25000.10.Cultivar_A-Cultivar_OutGroup-Wild_Mix_GroupA.txt \
                     10
                     
    Example: Rscript IntroPlot.R \
                     Dxy \
                     Dxy.50000.25000.10.Cultivar_A-Cultivar_OutGroup.txt \
                     10
                     
## 2) Example plot

High D-statistic value and low RND and Dxy values indicate candidate region for introgression events.

![그림3](https://user-images.githubusercontent.com/49300659/63839467-5180b400-c9ba-11e9-8468-196bcb3737f4.png)

                     
# Calculation Flow

Calculation flow of INTRO is as follows.

![그림1](https://user-images.githubusercontent.com/49300659/63830448-1fb22200-c9a7-11e9-86f5-ba709246719c.jpg)


# Requirement

Python 3.0 program and numpy library are requiered for calculation.

R program and ggplot library are required for generating result plots.

# Contact

jaeyoonkim72@gmail.com

# License

INTRO is registered with the Korean Copyright Commission under accession umber C-2017-024343.

Source code in INTRO is publicly available and anyone can modify and redistribute it.

# Reference

Ai, H., Fang, X., Yang, B., Huang, Z., Chen, H., Mao, L., ... & Yang, J. (2015). Adaptation and possible ancient interspecies introgression in pigs identified by whole-genome sequencing. Nature genetics, 47(3), 217.

Feder, J. L., Xie, X., Rull, J., Velez, S., Forbes, A., Leung, B., ... & Aluja, M. (2005). Mayr, Dobzhansky, and Bush and the complexities of sympatric speciation in Rhagoletis. Proceedings of the National Academy of Sciences, 102(suppl 1), 6573-6580.

Nei, M., & Li, W. H. (1979). Mathematical model for studying genetic variation in terms of restriction endonucleases. Proceedings of the National Academy of Sciences, 76(10), 5269-5273.

Patterson, N., Moorjani, P., Luo, Y., Mallick, S., Rohland, N., Zhan, Y., ... & Reich, D. (2012). Ancient admixture in human history. Genetics, 192(3), 1065-1093.

Rosenzweig, B. K., Pease, J. B., Besansky, N. J., & Hahn, M. W. (2016). Powerful methods for detecting introgressed regions from population genomic data. Molecular ecology, 25(11), 2387-2397.
