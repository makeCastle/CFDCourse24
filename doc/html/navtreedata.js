/*
 @licstart  The following is the entire license notice for the JavaScript code in this file.

 The MIT License (MIT)

 Copyright (C) 1997-2020 by Dimitri van Heesch

 Permission is hereby granted, free of charge, to any person obtaining a copy of this software
 and associated documentation files (the "Software"), to deal in the Software without restriction,
 including without limitation the rights to use, copy, modify, merge, publish, distribute,
 sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all copies or
 substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
 BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

 @licend  The above is the entire license notice for the JavaScript code in this file
*/
var NAVTREE =
[
  [ "cfd24", "index.html", [
    [ "Установка и сборка проекта", "install.html", [
      [ "Подготовка", "install.html#prep", null ],
      [ "VisualStudio", "install.html#vs_build", null ],
      [ "VSCode", "install.html#vscode_build", null ]
    ] ],
    [ "Описание проекта", "descr.html", null ],
    [ "Работа с инфраструктурой проекта", "howto.html", [
      [ "Git", "howto.html#git", [
        [ "Основные команды", "howto.html#autotoc_md0", null ],
        [ "Порядок работы с репозиторием CFDCourse", "howto.html#autotoc_md1", null ]
      ] ],
      [ "Сборка и запуск", "howto.html#autotoc_md2", [
        [ "Запуск конкретного теста", "howto.html#autotoc_md3", null ],
        [ "Сборка релизной версии", "howto.html#release_build", null ]
      ] ],
      [ "Paraview", "howto.html#paraview", [
        [ "Отображение одномерных графиков", "howto.html#paraview-1d", null ],
        [ "Отображение изолиний для двумерного поля", "howto.html#autotoc_md4", null ],
        [ "Отображение двумерного поля в 3D", "howto.html#paraview-2d", null ],
        [ "Отображение числовых данных для точек и ячеек", "howto.html#autotoc_md5", null ]
      ] ],
      [ "CMake", "howto.html#autotoc_md6", [
        [ "Добавление файла в проект", "howto.html#autotoc_md7", null ]
      ] ],
      [ "Построение одномерных графиков в логарифмических осях", "howto.html#autotoc_md8", [
        [ "Excel", "howto.html#autotoc_md9", null ],
        [ "LibreCalc", "howto.html#autotoc_md10", null ],
        [ "Octave/Matlab", "howto.html#autotoc_md11", null ],
        [ "Python", "howto.html#autotoc_md12", null ]
      ] ]
    ] ],
    [ "Краткое содержание лекций", "notes.html", [
      [ "1. Введение (02.09)", "notes.html#autotoc_md13", null ],
      [ "2. Уравнение Пуассона (09.09)", "notes.html#autotoc_md14", [
        [ "Постановка задачи", "notes.html#autotoc_md15", null ],
        [ "Метод решения", "notes.html#autotoc_md16", [
          [ "Нахождение численного решения", "notes.html#poisson1d_fdm", null ],
          [ "Практическое определения порядка аппроксимации", "notes.html#autotoc_md17", null ]
        ] ],
        [ "Программная реализация", "notes.html#test_poisson1", [
          [ "Функция верхнего уровня", "notes.html#autotoc_md18", null ],
          [ "Детали реализации", "notes.html#autotoc_md19", null ]
        ] ]
      ] ],
      [ "3. Двухслойные схемы для нестационарных уравнений (16.09)", "notes.html#autotoc_md20", [
        [ "Определение", "notes.html#autotoc_md21", [
          [ "Явная схема", "notes.html#autotoc_md22", null ],
          [ "Неявная схема", "notes.html#autotoc_md23", null ],
          [ "Схема Кранка–Николсон", "notes.html#autotoc_md24", null ],
          [ "Обобщённая двухслойная схема", "notes.html#autotoc_md25", null ]
        ] ],
        [ "Дискретизация по времени как итерационный процесс", "notes.html#autotoc_md26", [
          [ "Двухслойный итерационный процесс", "notes.html#autotoc_md27", null ],
          [ "Устойчивость итерационного процесса", "notes.html#ScalarIter", null ],
          [ "Источники возмущений", "notes.html#autotoc_md28", null ]
        ] ],
        [ "Методы исследования устойчивости расчётных схем", "notes.html#autotoc_md29", [
          [ "Матричный метод", "notes.html#autotoc_md30", [
            [ "Явная схема для нестационарного уравнения диффузии", "notes.html#NonstatExpDiff", null ],
            [ "Неявная схема для нестационарного уравнения диффузии", "notes.html#NonstatImpDiff", null ]
          ] ],
          [ "Метод дискретных возмущений", "notes.html#autotoc_md31", [
            [ "Явная схема против потока для уравнения переноса", "notes.html#NonstatExpTran", null ]
          ] ],
          [ "Метод Неймана", "notes.html#autotoc_md32", [
            [ "Неявная противопотоковая схема для уравнения переноса", "notes.html#NonstatImpConv", null ],
            [ "Противопотоковая схема Кранка-Николсон для уравнения переноса", "notes.html#NonstatCNConv", null ],
            [ "Явная схема для уравнения нестационарной конвекции-диффузии", "notes.html#NonstatExpConvDiff", null ],
            [ "Неявная схема для уравнения нестационарной конвекции-диффузии", "notes.html#NonstatImpConvDiff", null ]
          ] ],
          [ "Общие рекомендации к выбору устойчивых расчётных схем", "notes.html#autotoc_md33", null ]
        ] ],
        [ "Программная реализация схемы для уравнения переноса", "notes.html#autotoc_md34", [
          [ "Постановка задачи", "notes.html#autotoc_md35", null ],
          [ "Функция верхнего уровня", "notes.html#autotoc_md36", null ],
          [ "Расчётные функции", "notes.html#autotoc_md37", [
            [ "Явная схема", "notes.html#autotoc_md38", null ],
            [ "Неявная схема", "notes.html#autotoc_md39", null ],
            [ "Схема Кранка-Николсон", "notes.html#autotoc_md40", null ]
          ] ],
          [ "Анализ результатов работы", "notes.html#autotoc_md41", null ]
        ] ]
      ] ]
    ] ],
    [ "Задания для самостоятельной работы", "tasks.html", [
      [ "Лекция 1 (02.09)", "tasks.html#autotoc_md42", null ],
      [ "Лекция 2 (09.09)", "tasks.html#autotoc_md43", null ],
      [ "Лекция 3 (16.09)", "tasks.html#autotoc_md44", [
        [ "Постановка задачи", "tasks.html#autotoc_md45", [
          [ "Тестовый пример 1", "tasks.html#autotoc_md46", null ],
          [ "Тестовый пример 2", "tasks.html#autotoc_md47", null ]
        ] ],
        [ "Расчётная схема", "tasks.html#autotoc_md48", null ]
      ] ]
    ] ],
    [ "Namespaces", "namespaces.html", [
      [ "Namespace List", "namespaces.html", "namespaces_dup" ],
      [ "Namespace Members", "namespacemembers.html", [
        [ "All", "namespacemembers.html", null ],
        [ "Functions", "namespacemembers_func.html", null ]
      ] ]
    ] ],
    [ "Classes", "annotated.html", [
      [ "Class List", "annotated.html", "annotated_dup" ],
      [ "Class Index", "classes.html", null ],
      [ "Class Hierarchy", "hierarchy.html", "hierarchy" ],
      [ "Class Members", "functions.html", [
        [ "All", "functions.html", null ],
        [ "Functions", "functions_func.html", null ],
        [ "Typedefs", "functions_type.html", null ]
      ] ]
    ] ],
    [ "Files", "files.html", [
      [ "File List", "files.html", "files_dup" ],
      [ "File Members", "globals.html", [
        [ "All", "globals.html", null ],
        [ "Macros", "globals_defs.html", null ]
      ] ]
    ] ]
  ] ]
];

var NAVTREEINDEX =
[
"annotated.html",
"namespacemembers.html"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';