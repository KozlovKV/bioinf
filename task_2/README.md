# Козлов Кирилл Вячеславович - Биоинформатика-1 - задание 2
## Ген человека
Возьму уже использованный ранее ген [CRYGC](https://www.omim.org/entry/123680), кодирующий белок для хрусталика глаза.
- [NCBI](https://www.ncbi.nlm.nih.gov/gene/1420)
- [FASTA-файл](./mRNAs/crygc-hs.fna)

## Гомологи
Название на латыни | Название на русском | NCBI | Файл
:-- | :-- | :-: | :-: | 
Chinchilla lanigera | Длиннохвостая шиншила | [XM_005397113](https://www.ncbi.nlm.nih.gov/nuccore/XM_005397113) | [./mRNAs/1.fna](./mRNAs/1.fna)
Peromyscus maniculatus bairdii | Полевая оленья мышь | [XM_006974986.2](https://www.ncbi.nlm.nih.gov/nuccore/XM_006974986.2) | [./mRNAs/2.fna](./mRNAs/2.fna)
Lepus europaeus | Европейский заяц | [XM_062201647.1](https://www.ncbi.nlm.nih.gov/nuccore/XM_062201647.1) | [./mRNAs/3.fna](./mRNAs/3.fna)
Cynocephalus volans | Филипинский летающий лемур | [XM_063102653.1](https://www.ncbi.nlm.nih.gov/nuccore/XM_063102653.1) | [./mRNAs/4.fna](./mRNAs/4.fna)
Manis javanica | Индокитайский ящер (яванский ящер, яванский панголин) | [XM_036994993.2](https://www.ncbi.nlm.nih.gov/nuccore/XM_036994993.2) | [./mRNAs/5.fna](./mRNAs/5.fna)
Zalophus californianus | Калифорнийский морской лев | [XM_027589606.2](https://www.ncbi.nlm.nih.gov/nuccore/XM_027589606.2) | [./mRNAs/6.fna](./mRNAs/6.fna)
Vombatus ursinus | Короткошёрстый вомбат (эндемик астралии) | [XM_027873732](https://www.ncbi.nlm.nih.gov/nuccore/XM_027873732) | [./mRNAs/7.fna](./mRNAs/7.fna)
Anguilla rostrata | Американский угорь | [XM_064311543](https://www.ncbi.nlm.nih.gov/nuccore/XM_064311543) | [./mRNAs/8.fna](./mRNAs/8.fna)
Austrofundulus limnaeus | Нет собственного названия, но относятся к киллифишам | [XM_014020056](https://www.ncbi.nlm.nih.gov/nuccore/XM_014020056) | [./mRNAs/9.fna](./mRNAs/9.fna)
Alosa alosa | Аллис шад | [XM_048238469.1](https://www.ncbi.nlm.nih.gov/nuccore/XM_048238469.1) | [./mRNAs/10.fna](./mRNAs/10.fna)

(Все цепочки вместе для выравнивания собрал [тут](./mRNAs/mRNAs.fna))

## Множественное выравнивание
Использовал Clustal Omega

[Файл выравнивания](./mRNAs/alignment.aln-clustal_num)

Нуклеотидные цепочки получились самой разной длины у разных видов и были выравнены очень плохо, поэтому для получения хоть сколько-нибудь интересных результатов будем проводить выравнивания по белоквым цепочкам

---

В папке `./proteins` собрал данные по кодируемым белкам:
- [Белок человека](./proteins/hs.faa)
- [Все белки для выравнивания](./proteins/proteins.faa)
- [Выравнивание](./proteins/alignment.aln-clustal_num)

### Анализ
Здесь уже мы можем наблюдать достаточно большие консервативные участки, что позволяет нам выявить эволюционно важные участки данного белка, которые почти без изменений передаются уже очень долго.

При этом в последних нескольких цепочках есть ощутимые отличия от прочих, что логично, так как это рыбы, а не млекопитающие

## Общий таксон
Самый нижний общий таксон - *Euteleostomi (Костные позвоночные)* - включает в себя почти 90% всех позвоночных. Основные представители - рыбы и млекопитающие