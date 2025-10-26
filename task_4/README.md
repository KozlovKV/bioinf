# *Козлов Кирилл*, Биоинформатика-2, Задание 4 - Визуализация структуры белка
## Программа для визуализации
В качестве программы взял [Jmol](https://jmol.sourceforge.net/)

## Выбранный белок
В качестве белка выбрал белок гемоглобина [2PGH](https://www.rcsb.org/structure/2PGH)

## Выполнение работы
Открываем файл `jmol.bat` - у нас открывается приложение.

### Загрузка белка
В меню `File -> Get PDB` вводим код белка `2PGH`, после чего он автоматически подгрузится в программу

![](./screenshots/1_get.png)

### Визуализация
#### Wireframe
Нажимем `Display -> Bond -> Wireframe` либо ПКМ: `Style -> Scheme -> Wireframe`, после чего можем наблюдать структуру связей атомов в белке:

<img src="./screenshots/2_wireframe.png" width="45%" />
<img src="./screenshots/3_wireframe_res.png" width="45%" />

#### Backbone
ПКМ: `Style -> Scheme -> Trace`:

<img src="./screenshots/4_trace.png" width="45%" />
<img src="./screenshots/5_trace_res.png" width="45%" />

#### Spacefill
Атомы с полями

ПКМ: `Style -> Scheme -> CPK Spacefill`:

<img src="./screenshots/6_spacefill.png" width="45%" />
<img src="./screenshots/7_spacefill_res.png" width="45%" />

#### Ribbons
Ленточная вторичная структура

ПКМ: `Style -> Structures -> Ribbons`

<img src="./screenshots/8_ribbons.png" width="45%" />
<img src="./screenshots/9_ribbons_res.png" width="45%" />

#### Molecular surface
Позволит также отобразить поверхность белка. Наложу поверх ленточной структуры

ПКМ: `Surfaces -> Molecular surface`

<img src="./screenshots/10_surface.png" width="45%" />
<img src="./screenshots/11_surface_res.png" width="45%" />

### Окраска
#### CPK
На самом деле, это та модель окраски, которая видна по умолчанию, но если вдруг к ней надо вернуться после отображения других видов, то надо включить отображения атомов + ПКМ: `Color -> Atoms -> By Scheme -> Element (CPK)`, после чего сможем наблюадть окраску атомов (или их групп) в зависиомсоти от элементов:

<img src="./screenshots/12_color_CPK.png" width="45%" />
<img src="./screenshots/13_color_CPK_res.png" width="45%" />

#### Окраска по доменам
Включается через ПКМ: `Color -> By Scheme -> Chain`

<img src="./screenshots/14_doamin.png" width="45%" />

![](./screenshots/15_domain_res_1.png)

Можно отчётливо наблюдать 4 домена

## Снимки публикационного качества
С папке [`./export_images/`](./export_images/) можно посмотреть несколько изображений повышенного качества. Также приложу их здесь

![](./export_images/colored_by_chain_cartoon_rocket_with_dot_surface.jpg)

![](./export_images/colored_by_chain_ribbons_with_molecular_surface.jpg)