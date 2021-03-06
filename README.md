# InterAtomicInteraction
* Проект использует систему сборки _Cmake_ в режиме _out of source build_, т.е. перед запуском в корневой директории необходимо создать директорию `build/` и собрать проект под свою системную среду командой `cmake -G [среда разработки | Makefile] ..`  
Также в проекте используются сторонние библиотеки: `OpenMP` для параллельной разработки, `Catch2` для юнит-тестирования, `nlohmann::json` для парсинга `json`-файлов  

* Ломать зависимости в проекте рекомендуется ___только на свой страх и риск___

* `src` - директория с исполняющими файлами проекта  
* `test` - юнит тесты базовых структур проекта  

* `input_data.json` – входной файл с табличными параметрами для функционала минимизации и начальными приближениями параметров потенциалов для процедуры оптимизации. Табличные параметры заданы в строгом порядке _a_, _Ecoh_, _B_, _C_11_, _C_12_, _C_44_, _E_sol_, _E_in_dim_, _E_on_dim_, _E_coh_A_dim_ как и начальные приближения для параметров потенциалов _A_0_, _A_1_, _ksi_, _p_, _q_, _r_0_  
* `vector3D.hpp` - класс трехмерного вектора  
* `matrix3D.hpp` - класс трехмерной матрицы деформации  
* `atom.h, atom.cpp` - класс атома в трехмерном пространстве  
* `lattice.h, lattice.cpp` - класс решетки ГЦК  
* `table_estimator.h, table_estimator.cpp` - решатель прямого хода задачи: подбор табличных параметров по известным параметрам потенциала (известным - в смысле рассматриваемым на текущей итерации процесса оптимизации в обратном ходе задачи)  
* `optimizer.h, optimizer.cpp` - решатель обратного хода задачи: оптимизация функции ошибки для табличных параметров по параметрам потенциала (ошибка между посчитанными параметрами в прямом ходе задачи и истинными значениями из таблицы)  
