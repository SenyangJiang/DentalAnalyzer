# Developer's Guide
## OS
We will develop this program under Windows 10.
## Libraries
Below lists the libraries that I used:
- [CGAL](https://www.cgal.org/): Version 5.1.1
- [boost](https://www.boost.org/): Dependency of CGAL, Version 1.75
- [QT](https://www.qt.io/cs/c/?cta_guid=074ddad0-fdef-4e53-8aa8-5e8a876d6ab4&signature=AAH58kGRjlak3s62ZJ7smf4xZ0fAiA2iow&pageId=12602948080&placement_guid=99d9dd4f-5681-48d2-b096-470725510d34&click=06f85479-4897-4374-a7e5-f0c3ef28760b&hsutk=7197fbf4fef3f20fa6a3170fadcede60&canon=https%3A%2F%2Fwww.qt.io%2Fdownload-open-source&portal_id=149513&redirect_url=APefjpFZbAR6UAJKFS84MfBMZQCPe_fllUoYe3yFgJpndM5YNWG5-cguvwbHULwltllTpsKxphwf7yPWs2cFA5l6G6yW-Oulq5U-u-B-OAuuDKL0skaVHT6KUT6RQyf5nO4ltuPnL-kX&__hstc=152220518.7197fbf4fef3f20fa6a3170fadcede60.1665196588673.1665196588673.1665196588673.1&__hssc=152220518.1.1665196588673&__hsfp=2042877427&contentType=standard-page): IDE for development, convenient to create UI, Version 5.15.2 (msvc2019_64)

I will attach the CGAL and boost library here. Download and put them in `C:/dev`

https://drive.google.com/file/d/1dbSDf-wNOQJwgzVqnwHMa4O_Jo3yn6M1/view?usp=sharing

You can add the QT library in the IDE. Make sure it is QT5.15 msvc2019_64.

## Files
**A general overview:**
- `CMakelists.txt`: information about library linking in compilation
- `main.cpp`: generated files for UI where I haven't touched
- `mainwindow.ui`: If you double click this file, QT takes you to the UI page where you can modify the UI
- `mainwindow.h`, `mainwindow.cpp`: generated files for UI. The UI logics are defined here. For example, if you add a button, the button logic are defined in these files. To learn more about UI logic, Google "QT slots and signals".
- `Analyzer.cpp`, `Analyzer.h`: These two files defines Analyzing functions in the Analyzer class.
- `parameter.h`: parameter to the Analyzer object
- `objects.h`: defines some CGAL objects
- `result.h`, `result.cpp`: class definitions for storing Analyzer results.
- `utility.h`, `utility.cpp`: some utility functions

## Example workflow and how it's related to the code
1. User opens `DentalAnalyzer.exe`, they will see the UI, which is defined in `mainwindow.ui`.
2. User click buttons to input their 3D model files, the buttons' logic are defined in `mainwindow.h` and `mainwindow.cpp`. Most of these buttons will read the input string and store it in a parameter object. The parameter class is defined in `parameter.h`.
3. When the user click "Analyze", it will trigger the `on_pushButtonAnalyze_clicked()` function in `mainwindow.cpp`. This function will create a analyer object and run its `analyzer()` function in a **new** thread. If you directly run it, the mainwindow will freeze because the main thread is busy running the analyze function.
4. The `Analyzer::analyze()` function will read in the models, and run the compute functions. When running the functions, the result will be saved to the `student_result` variable in Analyzer object.
5. When you export the result as csv using the button, the result object from the previous step will be used and write data to a csv file.