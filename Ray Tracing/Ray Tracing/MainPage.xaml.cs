using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices.WindowsRuntime;
using Windows.Foundation;
using Windows.Foundation.Collections;
using Windows.UI.Xaml;
using Windows.UI.Xaml.Controls;
using Windows.UI.Xaml.Controls.Primitives;
using Windows.UI.Xaml.Data;
using Windows.UI.Xaml.Input;
using Windows.UI.Xaml.Media;
using Windows.UI.Xaml.Navigation;

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
        public MainPage()
        {
            this.InitializeComponent();
            timer = new DispatcherTimer() { Interval = new TimeSpan(0, 0, 1) };
            timer.Tick += Timer;
            timer.Start();
        }

        private void ToogleSettings(object sender, RoutedEventArgs e)
        {
            SettingsMenu.IsPaneOpen = !SettingsMenu.IsPaneOpen;
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
        }
    }
}
