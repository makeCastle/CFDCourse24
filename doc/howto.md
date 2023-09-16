@page howto Работа с инфраструктурой проекта

# Git {#git}

## Основные команды

Все команды выполнять в теримнале (git bash для виндоус),
находясь в корневой папке проета CFDCourse24.

- Для **смены директории** использовать команду cd. Например, находясь в папке A перейти в папку A/B/C
  @code
  > cd B/C
  @endcode

- **Подняться** на директорию выше
  @code
  > cd ..
  @endcode

- **Просмотр статуса** текущего репозитория: текущую ветку, все изменённые файлы и т.п.
  @code
  > git status
  @endcode

- **Сохранить и скоммитить** изменения в текущую ветку
  @code
  > git add .
  > git commit -m "message"
  @endcode

  "message" -- произвольная информация о текущем коммите, которая будет приписана к этому коммиту

- **Переключиться на ветку** main
  @code
  > git checkout main
  @endcode

  работает только в том случае, если все файлы скоммичены и статус ветки 'Up to date'

- **Создать новую ветку** ответвлённую от последнего коммита текущей ветки и переключиться на неё
  @code
  > git checkout -b new-branch-name
  @endcode

  new-branch-name -- имя новой ветки. Пробелы не допускаются

  Эта комманда работает даже если есть нескоммиченные изменения. 
  Если необходимо скоммитить изменеия в новую ветку, сразу за этой командой нужно вызвать
  @code
  > git add .
  > git commit -m "message"
  @endcode

- **Сбросить** все нескоммиченные изменения. Вернуть файлы в состояние последнего коммита
  @code
  > git reset --hard
  @endcode

  Все изменения будут утеряны

- **Получить последние изменения** из удалённого хранилища с обновлением текущей ветки
  @code
  > git pull
  @endcode
  Работает только если статус текущей ветки 'Up to date'.\n
  Если требуется получить изменения, но не обновлять локальную ветку:
  @code
  > git fetch
  @endcode
  Обновленная ветка будет доступна по имени origin/{имя ветки}.

- **Просмотр истории** коммитов в текущей ветке (последний коммит будет наверху)
  @code
  > git log
  @endcode

- **Просмотр** актуального состояние дерева репозитория в gui режиме
  @code
  > git gui
  @endcode
  Далее в меню Repository->Visualize all branch history.
  В этом же окне можно посмотреть изменения файлов по сравнению с последним коммитом.

  Альтернативно, при работе в виндоус можно установить программу GitExtensions и работать в ней.
  
## Порядок работы с репозиторием CFDCourse

Основная ветка проекта -- main. После каждой лекции (в течении 1-2 дней) в эту ветку будет отправлен коммит с сообщением "after-lect{index}".
Этот коммит будет содержать краткое \ref notes "содержание лекции",
\ref tasks "задание" по итогам лекции и необходимые для этого задания изменения кода.

Перед лекцией в эту ветку будет отправлен коммит с сообщением "before-lect{index}".
Этот коммит содержит изменения кода для работы на лекции.

Таким образом, **после лекции** необходимо выполнить следующие команды (находясь в ветке main)
@code
> git reset --hard  # очистить локальную копию от изменений, сделанных на лекции (если они не представляют ценности)
> git pull
@endcode

**Перед началом лекции**, если была сделана какая то работа по заданиям
@code
> git checkout -b work-lect{index}    # создать локальную ветку, содержащую задание
> git add .
> git commit -m "{свой комментарий}"  # скоммитить свои изменения в эту ветку
> git checkout main                   # вернуться на ветку main
> git pull                            # получить изменения
@endcode

Даже если задание выполнено не до конца, вы в любой момент можете переключиться на ветку с заданием и его доделать
@code
> git checkout work-lect{index}
@endcode

Если ничего не было сделано (или все изменения не представляют ценности), можно повторить алгоритм "после лекции".

# Сборка и запуск

## Запуск конкретного теста

По умолчанию программа cfd_test прогоняет все объявленные в проекте тесты. Иногда может возникнуть необходимость
запустить только конкретный тест в целях отладки или проверки.
Для этого нужно передать программе аргумент с тегом для этого теста.

Тег для теста -- это второй аргумент в макросе \c TEST_CASE, записанный в квадратных скобках.
Добавлять нужно вместе со скобками. Например, "[ping]".

Чтобы добавить аргумент в VisualStudio, необходимо в контекстном меню проекта cfd_test выбрать опции отладки
\image html win_debug_args_1.png
и там в поле Аргументы прописать нужный тэг.
\image html win_debug_args_2.png

В VsCode аргументы нужно добавлять в файле \c .vscode/launch.json в поле args в кавычках (см. картинку с \ref vscode_build "настройками launch.json").

## Сборка релизной версии {#release_build}

TODO

# Paraview {#paraview}

## Отображение одномерных графиков {#paraview-1d}

Заданные на сетке данные паравью показывает цветом.
Поэтому при загрузке одномерных графиков можно видеть следующую картинку

\image html howto_paraview_1d_1.png width=600

В первую очередь нужно развернуть изображение в плоскость xy
\image html howto_paraview_1d_2.png width=400

Далее, для того, что бы данные отображались в качестве значения по оси ординат, к загруженному файлу необходимо
    1. применить фильтр WarpByScalar (В меню Filters->Alphabetical->Warp By Scalar)
    2. в меню настройки фильтра указать поле данных, для отображения (numerical в примере ниже)
    3. И настроить нормаль, вдоль которой будут проецироваться данные (в нашем случае ось y)
\image html howto_paraview_1d_3.png width=600

Чтобы настроить цвет и толщину линии нужно
    1. Включить подробные опции фильтра
    2. Сменить стиль на Solid Color
    3. В меню Edit выбрать желаемый цвет
    4. В строке Line Width указать толщину линии
\image html howto_paraview_1d_4.png width=600

Настрока масштабов и отображение осей координат:
    1. Отметье подробные настройки фильтра
    2. В поле Transforming/Scale Установите желаемые масштабы (в нашем случае растянуть в два раза по оси x)
    3. Установите галку на отображение осей
    4. откройте меню натройки осей
    5. В нём включите подробные настроки
    6. И также поставьте растяжение осей

В случае, если мастабировать график не нужно, достаточно выполнить шаг 3.
\image html howto_paraview_1d_5.png width=1000

Если требуется нарисовать рядом несколько графиков для разных данных из одного файла,
примените фильтр Warp By Scalar для этого файла ещё раз, изменив поле Scalars в настройке фильтра.
Для наглядности измените имя узла в Pipeline Browser на осмысленные
\image html howto_paraview_1d_6.png width=800

В случае, если исходный файл был изменён, нужно в контекстном меню узла соответствующего файла
выбрать Reload Files (или нажать F5). Если те же самые фильтры нужно применить для просмотра другого файла
нужно в этом меню нажать Change File.
\image html howto_paraview_1d_7.png width=300

## Отображение изолиний для двумерного поля

1. Нажмите иконку Contour (или Filters/Contour) (1 на рисунке)
   В настройках фильтра Contour by выберитее данные, по которым нужно строить изолинии.
2. В настройках фильтра удалите все существующие записи о значениях для изолиний (2 на рисунке)
3. Добавьте равномерные значения (3 на рисунке). В появившемся меню установите необходимое количество изолиний и их диапазон.
4. Если необходимо, включите одновременное отображения цветного поля и изолиний (4 на рисунке).

\image html howto_paraview_isolines_1.png width=800

В случае, если нужно сделать изолинии одного цвета, установите поле Coloring/Solid color в 
настройках фильтра. Там же в меню Edit можно выбрать цвет.

Для установления толщины линии включите подробные настройки и найдите там опцию Styling/Line Width.

\image html howto_paraview_isolines_2.png width=800

## Отображение двумерного поля в 3D
По аналогии с \ref paraview-1d "одномерным графиком", двумерные поля так же
можно отобразить, проектируя данные на геометрическую координату для получения
объёмного графика.

1. Включите фильтр Filters/Warp By Scalar
2. В настройках фильтра установите данные, которые будут проектироваться на координату z
3. Установите нормаль для проецирования (ось z)
4. Если нужно, выберите масштабирования для этой координаты
5. После нажатия Apply включите трёхмерное отображение
6. Если данные не видно, обносите экран.

\image html howto_paraview_2d_as_3d.png width=800

## Отображение числовых данных для точек и ячеек
1. Включить режим выделения точек или ячеек (иконка (1 на рисунке) или горячие клавиши s, d)\n
2. Выделить мышкой интересующую область
3. В окне Find data (или Selection Inspector для старых версий Paraview) отметить поле, которое должно отображаться 
   в центрах ячеек и в точках (2 на рисунке). Если такого окна нет, включить его из основного меню View.

\image html howto_paraview_show_labels.png width=800


# CMake
## Добавление файла в проект
Для добавления нового ресурсного файла наобходимо прописать его в соответствующий \c CMakeLists.txt.
- Если файл добавляется в библиотеку \c cfd, то \c src/cfd24/CMakeLists.txt
- Если файл добавляется в тестовое приложение, то \c src/test/CMakeLists.txt

Заголовочные файлы (расширение hpp) прописывать в список \c HEADERS,
компилируемые (расширение cpp) прописывать в \c SRC.

При работе с VisualStudio (если он запущен) после изменения cmake-файла необходимо построить проект \c ZERO_CHECK
вызвав контекстное меню и вызвав \c Build для обновления файлового дерева.
\image html win_zero_check_build.png width=400
В результате должен появится диалог с
предложение обновить проект. Надо нажать Reload All (Обновить всё).