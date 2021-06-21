using System;
using System.Drawing;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices.WindowsRuntime;
using Windows.Foundation;
using Windows.Foundation.Collections;
using Windows.Storage.Pickers;
using Windows.UI;
using Windows.UI.Input;
using Windows.UI.Popups;
using Windows.UI.Xaml;
using Windows.UI.Xaml.Controls;
using Windows.UI.Xaml.Controls.Primitives;
using Windows.UI.Xaml.Data;
using Windows.UI.Xaml.Input;
using Windows.UI.Xaml.Media;
using Windows.UI.Xaml.Media.Imaging;
using Windows.UI.Xaml.Navigation;
using Image = Windows.UI.Xaml.Controls.Image;
using Page = Windows.UI.Xaml.Controls.Page;
using Point = Windows.Foundation.Point;
using Windows.UI.Xaml.Shapes;
// Документацию по шаблону элемента "Пустая страница" см. по адресу https://go.microsoft.com/fwlink/?LinkId=402352&clcid=0x419

namespace Ray_Tracing
{
    /// <summary>
    /// Пустая страница, которую можно использовать саму по себе или для перехода внутри фрейма.
    /// </summary>
    public sealed partial class MainPage : Page
    {
        DispatcherTimer timer;
        int time = 0;
        List<Circle> list = new List<Circle>();
        int M = 0;		//x column
        int N = 0;		//y string
        int L = 0;		//Монте-Карло (кол-во траекторий)
        int scatt = 0;	// акты рассеяния
        int c = 0;		//скорость звука
        int cancelClick = 0, historyClick = 0, sizeClick = 0, timeResolve = 0;
        int sigmatimes = 0;	//множитель
        double Mu = 0;	//коэф. затухания
        double h = 0.4;	//шаг
        double PI = 3.1415; //пи
        double maximum = 0; //максимум

        public MainPage()
        {
            this.InitializeComponent();
            timer = new DispatcherTimer() { Interval = new TimeSpan(0, 0, 1) };
            timer.Tick += Timer;
            timer.Start();
            Windows.UI.ViewManagement.ApplicationView appView = Windows.UI.ViewManagement.ApplicationView.GetForCurrentView();
            appView.Title = "Моделирование нестационарных процессов зондирования неоднородных рассеивающих сред";
        }

        [DllImport("RayTracindDLL.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern int main();

        private void ToogleSettings(object sender, RoutedEventArgs e)
        {
            SettingsMenu.IsPaneOpen = !SettingsMenu.IsPaneOpen;
        }

        private void ToogleParameters(object sender, RoutedEventArgs e)
        {
            string text = "Длина графика: " + M.ToString() + ";\n" +
                "Ширина графика: " + N.ToString() + ";\n" +
                "Монте-Карло (кол-во траекторий): " + L.ToString() + ";\n" +
                "Акты рассеяния: " + scatt.ToString() + ";\n" + 
                "Скорость звука: " + c.ToString() + ";\n" + 
                "Множитель: " + sigmatimes.ToString() + ";\n" +
                "Коэффициент затухания: " + Mu.ToString() + ";\n" +
                "Шаг: " + h.ToString() + ";\n" +
                "Пи: " + PI.ToString() + ";\n" +
                "Максимум: " + maximum.ToString() + ";\n";
            ContentDialog parameters = new ContentDialog()
            {
                Title = "Текущие параметры",
                Content = text,
                PrimaryButtonText = "ОК"
            };
            parameters.ShowAsync();
        }

        private void Timer(object sender, object e)
        {
            time++;
            string hour = (time / 3600).ToString();
            if (hour.Length < 2)
            {
                hour = "0" + hour;
            }
            string minute = (time / 60).ToString();
            if (minute.Length < 2)
            {
                minute = "0" + minute;
            }
            string second = (time % 60).ToString();
            if (second.Length < 2)
            {
                second = "0" + second;
            }
            Time.Text = "Текущая сессия: " + hour + ":" + minute + ":" + second;
            if (Settings.IsEnabled == false)
            {
                timeResolve += 1;
            }
            if (timeResolve == 60)
            {
                Settings.IsEnabled = true;
                timeResolve = 0;
                ProgressRing.IsActive = false;
                ProgressScope.Visibility = Visibility.Collapsed;
                string text = "Готово! Можете посмотреть результат обработки в истории";
                ContentDialog parameters = new ContentDialog()
                {
                    Title = "Завершение обработки",
                    Content = text,
                    PrimaryButtonText = "ОК"
                };
                parameters.ShowAsync();
            }
        }

        private void Grid_PointerMoved(object sender, PointerRoutedEventArgs e)
        {
            if (sizeClick != 0)
            {
                PointerPoint pp = e.GetCurrentPoint(relativeTo: Graph);
                Point point = pp.Position;
                int correct = Convert.ToInt32(xValue.Text) / 2;
                int xCoord = Convert.ToInt32(point.X) - correct;
                int yCoord = Convert.ToInt32(point.Y);
                Coordinates.Text = "Текущие координаты: {" + xCoord + " " + yCoord + "}";
            }
            else
            {
                Coordinates.Text = "Текущие координаты: ";
            }
        }

        private void setSpace_Click(object sender, RoutedEventArgs e)
        {
            Graph.Width = Convert.ToInt32(xValue.Text);
            Graph.Height = Convert.ToInt32(yValue.Text);
            setParam.IsEnabled = true;
            sizeClick = 1;
        }

        private void ValueTextChanging(TextBox sender, TextBoxTextChangingEventArgs args)
        {
            if (xValueCircle.Text == "")
            {
                xValueCircle.Text = "0";
            }
            else if (xValueCircle.Text == "-") 
            {
                xValueCircle.Text = "0";
            }
            else if (Convert.ToInt32(xValueCircle.Text) < 0)
            {
                return;
            }
            int pos = sender.SelectionStart;
            sender.Text = new String(sender.Text.Where(char.IsDigit).ToArray());
            sender.SelectionStart = pos;
        }

        private void ValueLostFocus(object sender, RoutedEventArgs e)
        {
            if (hValue.Text == "")
            {
                hValue.Text = "0,1";
            }
            if (sigmaValue.Text == "")
            {
                sigmaValue.Text = "1";
            }
            Scatter.Text = "Количество актов рассеяния: " + SValue.Text;
            if (xValue.Text == "")
            {
                xValue.Text = "400";
            }
            if (Convert.ToInt32(xValue.Text) > 1400)
            {
                xValue.Text = "1400";
            }
            if (Convert.ToInt32(xValue.Text) == 0)
            {
                xValue.Text = "400";
            }
            if (Convert.ToInt32(xValue.Text) % 100 != 0)
            {
                int value = Convert.ToInt32(xValue.Text) % 100;
                if (value >= 50)
                {
                    xValue.Text = (Convert.ToInt32(xValue.Text) + 100 - value).ToString();
                }
                else
                {
                    xValue.Text = (Convert.ToInt32(xValue.Text) - value).ToString();
                }
            }
            if (yValue.Text == "")
            {
                yValue.Text = "200";
            }
            if (Convert.ToInt32(yValue.Text) > 600)
            {
                yValue.Text = "600";
            }
            if (Convert.ToInt32(yValue.Text) == 0)
            {
                yValue.Text = "200";
            }
            if (Convert.ToInt32(yValue.Text) % 100 != 0)
            {
                int value = Convert.ToInt32(yValue.Text) % 100;
                if (value >= 50)
                {
                    yValue.Text = (Convert.ToInt32(yValue.Text) + 100 - value).ToString();
                }
                else
                {
                    yValue.Text = (Convert.ToInt32(yValue.Text) - value).ToString();
                }
            }
            if (LValue.Text == "")
            {
                LValue.Text = "1";
            }
            if (Convert.ToInt32(LValue.Text) > 1000)
            {
                LValue.Text = "1000";
            }
            if (Convert.ToInt32(LValue.Text) == 0)
            {
                LValue.Text = "1";
            }
            if (SValue.Text == "")
            {
                SValue.Text = "1";
                Scatter.Text = "Количество актов рассеяния: 1";
            }
            if (SValue.Text == "0")
            {
                SValue.Text = "1";
                Scatter.Text = "Количество актов рассеяния: 1";
            }
            if (Convert.ToInt32(SValue.Text) > 15)
            {
                SValue.Text = "15";
                Scatter.Text = "Количество актов рассеяния: 15";
            }
            if (sigmaValue.Text == "0")
            {
                sigmaValue.Text = "1";
            }
            if (Convert.ToInt32(sigmaValue.Text) > 5)
            {
                sigmaValue.Text = "5";
            }
            if (xValueCircle.Text == "")
            {
                xValueCircle.Text = "0";
            }
            if (Convert.ToInt32(xValueCircle.Text) > Convert.ToInt32(xValue.Text) / 2)
            {
                xValueCircle.Text = (Convert.ToInt32(xValue.Text) / 2).ToString();
            }
            else if (Convert.ToInt32(xValueCircle.Text) < -1 * (Convert.ToInt32(xValue.Text) / 2))
            {
                xValueCircle.Text = (-1 * Convert.ToInt32(xValue.Text) / 2).ToString();
            }
            if (yValueCircle.Text == "")
            {
                yValueCircle.Text = "0";
            }
            if (Convert.ToInt32(yValueCircle.Text) > Convert.ToInt32(yValue.Text))
            {
                yValueCircle.Text = yValue.Text;
            }
        }

        private void setParam_Click(object sender, RoutedEventArgs e)
        {
            M = Convert.ToInt32(xValue.Text);
            N = Convert.ToInt32(yValue.Text);
            L = Convert.ToInt32(LValue.Text);
            scatt = Convert.ToInt32(SValue.Text);
            sigmatimes = Convert.ToInt32(sigmaValue.Text);
            Mu = 0.018 * sigmatimes;
            h = Convert.ToDouble(hValue.SelectedValue);
            if (freshWater.IsChecked == true)
            {
                c = 1403;
            }
            else if (saltWater.IsChecked == true)
            {
                c = 1500;
            }
            else
            {
                c = 331;
            }
            maximum = 0.018 * sigmatimes;
            CurParam.IsEnabled = true;
            Save.IsEnabled = true;
        }

        private void saltWater_Click(object sender, RoutedEventArgs e)
        {
            freshWater.IsChecked = false;
            Air.IsChecked = false;
        }

        private void Air_Click(object sender, RoutedEventArgs e)
        {
            freshWater.IsChecked = false;
            saltWater.IsChecked = false;
        }

        private void freshWater_Click(object sender, RoutedEventArgs e)
        {
            saltWater.IsChecked = false;
            Air.IsChecked = false;
        }

        private void Draw_Click(object sender, RoutedEventArgs e)
        {
            int x = Convert.ToInt32(xValueCircle.Text) + (Convert.ToInt32(xValue.Text) / 2);
            int y = Convert.ToInt32(yValueCircle.Text);
            int d = Convert.ToInt32(dValueCircle.Text);
            int x_list = Convert.ToInt32(xValueCircle.Text) + (d / 2);
            int y_list = y + (d / 2);
            list.Add(new Circle(x_list, y_list)); //добавляем круг в список
            var ellipse = new Ellipse();
            ellipse.Fill = new SolidColorBrush(Windows.UI.Colors.White);
            ellipse.Width = d;
            ellipse.Height = d;
            Thickness thickness = new Thickness();
            thickness.Left = x;
            //thickness.Right = x;
            thickness.Top = y;
            //thickness.Bottom = y;
            ellipse.Margin = thickness;
            Graph.Children.Add(ellipse);
            Cancel.IsEnabled = true;
            cancelClick += 1;
        }

        private async void History_Click(object sender, RoutedEventArgs e)
        {
            Settings.IsEnabled = true;
            ProgressRing.IsActive = false;
            ProgressScope.Visibility = Visibility.Collapsed;
            if (historyClick == 1)
            {
                Graphics.Children.RemoveAt(2);
                historyClick -= 1;
            }
            historyClick += 1;
            FileOpenPicker openPicker = new FileOpenPicker();
            openPicker.ViewMode = PickerViewMode.Thumbnail;
            openPicker.SuggestedStartLocation = PickerLocationId.Desktop;
            openPicker.CommitButtonText = "Открыть";
            openPicker.FileTypeFilter.Add(".bmp");
            var file = await openPicker.PickSingleFileAsync();
            if (file == null)
            {
                historyClick -= 1;
                return;
            }
            Image image = new Image();
            Windows.UI.Xaml.Media.Imaging.BitmapImage bitmap = new Windows.UI.Xaml.Media.Imaging.BitmapImage();
            string uri = "ms-appx:///Assets/Image/" + file.DisplayName + ".bmp";
            bitmap.UriSource = new Uri(uri);
            image.Source = bitmap;
            image.Stretch = Stretch.UniformToFill;
            image.PointerMoved += Image_PointerMoved;
            Graphics.Children.Add(image);
        }

        private void Image_PointerMoved(object sender, PointerRoutedEventArgs e)
        {
            PointerPoint pp = e.GetCurrentPoint(relativeTo: Graph);
            Point point = pp.Position;
            int xCoord = Convert.ToInt32(point.X);
            int yCoord = Convert.ToInt32(point.Y);
            Coordinates.Text = "Текущие координаты: {" + xCoord + " " + yCoord + "}";
        }

        private void minusDegree_Click(object sender, RoutedEventArgs e)
        {
            xValueCircle.Text = (Convert.ToInt32(xValueCircle.Text) * -1).ToString();
        }

        private void Progress_Click(object sender, RoutedEventArgs e)
        {
            ProgressRing.IsActive = true;
            ProgressScope.Visibility = Visibility.Visible;
            //main();
            Settings.IsEnabled = false;
        }

        private void Cancel_Click(object sender, RoutedEventArgs e)
        {
            Graph.Children.RemoveAt(cancelClick - 1);
            cancelClick -= 1;
            if (cancelClick == 0)
            {
                Cancel.IsEnabled = false;
            }
        }
    }

    internal class Circle
    {
        public int X { get; set; }
        public int Y { get; set; }

        public Circle(int x, int y)
        {
            this.X = x;
            this.Y = y;
        }
    }
}
