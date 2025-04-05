package main

import (
	"fmt"
	"math"
)

// Функция
// 7x^2+3y^2+0.5xy-3x-5y+2
func f(x, y float64) float64 {
	return 7*x*x + 3*y*y + 0.5*x*y - 3*x - 5*y + 2
}

// Частная производная по x
func dfdx(x, y float64) float64 {
	return 14*x + 0.5*y - 3
}

// Частная производная по y
func dfdy(x, y float64) float64 {
	return 6*y + 0.5*x - 5
}

func findX(y float64) float64 {
	// выражаем x из ч.п. по y
	return (-0.5*y + 3) / 14
}

func findY(x float64) float64 {
	// выражаем y из ч.п. по x
	return (-0.5*x + 5) / 6
}

// Метод покоординатного спуска
// (x0, y0) - начальное приближение, tolerance - точность, maxIter - чтобы без бредика
func coordinateDescent(x0, y0, tolerance float64, maxIter int) (float64, float64) {
	x, y := x0, y0
	for iter := 0; iter < maxIter; iter++ {
		fmt.Println("Итерация ", iter, ":")
		// Минимизация по x
		newX := findX(y)

		fmt.Printf("x%d = %.4f\n", iter, newX)
		fmt.Printf("Текущая точка минимума M(%.4f, %.4f)\n", newX, y)

		if iter >= 1 {
			fmt.Printf("|f_x(M%d) - f_x(M%d)| = %.4f\n", iter, iter-1, f(newX, y)-f(x, y))
		}

		if math.Abs(f(newX, y)-f(x, y)) < tolerance {
			x = newX
			break
		}
		x = newX

		// Минимизация по y
		newY := findY(x)

		fmt.Printf("y%d = %.4f\n", iter, newY)
		fmt.Printf("Текущая точка минимума M(%.4f, %.4f)\n", x, newY)

		if iter >= 1 {
			fmt.Printf("|f_y(M%d) - f_y(M%d)| = %.4f\n", iter, iter-1, f(x, newY)-f(x, y))
		}

		if math.Abs(f(x, newY)-f(x, y)) < tolerance {
			y = newY
			break
		}
		y = newY
		fmt.Println("=================================================================================================")
	}
	fmt.Println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	return x, y
}

// Метод градиентного спуска
func gradDescent(x0, y0, alpha, tolerance float64, maxIter int) (float64, float64) {
	x, y := x0, y0
	for iter := 0; iter < maxIter; iter++ {
		fmt.Println("Итерация ", iter, ":")

		// Шаг против градиента по x
		newX := x - alpha*dfdx(x, y)
		// Шаг против градиента по y
		newY := y - alpha*dfdy(x, y)

		fmt.Printf("M%d_x Компонента Градиента = %.4f\n", iter, newX)
		fmt.Printf("M%d_y Компонента Градиента = %.4f\n", iter, newY)

		if iter >= 1 {
			fmt.Printf("|f(M%d) - f(M%d)| = %.4f\n", iter, iter-1, f(newX, newY)-f(x, y))
		}

		if math.Abs(f(x, y)-f(newX, newY)) < tolerance {
			break
		}
		x = newX
		y = newY

		fmt.Printf("Текущая точка минимума M(%.4f, %.4f)\n", x, y)
		fmt.Println("=================================================================================================")
	}
	fmt.Println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	return x, y
}

func calcRatios(x, y, dx, dy float64) []float64 {
	res := make([]float64, 0)
	// кэф перед x^2
	res = append(res, 7*dx*dx+3*dy*dy+0.5*dx*dy)
	// кэф перед x
	res = append(res, 7*(x*dx*2)+3*(y*dy*2)+0.5*(dx*y+x*dy)-3*dx-5*dy)
	// свободный кэф
	res = append(res, 7*x*x+3*y*y+0.5*(x*y)-3*x-5*y+2)

	return res
}

func findStep(ratios []float64) float64 {
	stepEqRatios := make([]float64, 0)
	// производная функции от шага
	stepEqRatios = append(stepEqRatios, 2*ratios[0])
	stepEqRatios = append(stepEqRatios, ratios[1])
	fmt.Printf("Вторая производная от функции шага: %.10f*h + %.10f\n", stepEqRatios[0], stepEqRatios[1])

	// приравняли к нулю и выразили шаг
	step := -stepEqRatios[1] / stepEqRatios[0]
	return step
}

func fastestDescent(x0, y0, tolerance float64, maxIter int) (float64, float64) {
	x, y := x0, y0
	ratios := make([]float64, 3)
	step := 0.0
	for iter := 0; iter < maxIter; iter++ {
		fmt.Println("Итерация ", iter, ":")
		// ч. п.
		derivX := dfdx(x, y)
		derivY := dfdy(x, y)
		fmt.Printf("Градиент grad = (%.6f*i, %.6f*j)\n", derivX, derivY)
		fmt.Printf("Модуль Градиента |grad| = %.6f\n", math.Abs(math.Sqrt(derivX*derivX+derivY*derivY)))

		if math.Abs(math.Sqrt(derivX*derivX+derivY*derivY)) < tolerance {
			break
		}
		// подставляем в исходную функцию уравнения для шага
		// x_new = x_prev - step * df/dx
		// y_new = y_prev - step * df/dy
		ratios = calcRatios(x, y, -derivX, -derivY)
		fmt.Printf("Подстановка в исходную функцию %.10f*h^2 + %.10f*h + %.10f\n", ratios[0], ratios[1], ratios[2])
		// считаем шаг через вторую производную + приравнивание к нулю
		step = findStep(ratios)
		fmt.Printf("Шаг h = %.4f\n", step)

		// новое приближение
		newX := x - step*derivX
		newY := y - step*derivY

		x = newX
		y = newY

		fmt.Printf("Текущая точка минимума M(%.4f, %.4f)\n", x, y)
		fmt.Println("=================================================================================================")
	}
	fmt.Println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	return x, y
}

func main() {
	x0, y0 := 2.0, -2.0 // Начальные приближения
	tolerance := 0.0001 // Точность
	maxIter := 1000     // Максимальное число итераций

	fmt.Println("Метод Покоординатного Спуска")
	fmt.Println("-------------------------------------------------------------------------------------------------")
	xMin, yMin := coordinateDescent(x0, y0, tolerance, maxIter)
	fmt.Printf("Минимум функции достигается в точке (x = %.4f, y = %.4f) со значением f(x,y) = %.4f\n", xMin, yMin, f(xMin, yMin))
	fmt.Println("=================================================================================================\n")

	fmt.Println("Метод Градиентного Спуска")
	fmt.Println("--------------------------------------------------------------------------------------------------")
	xMin, yMin = gradDescent(x0, y0, 0.1, tolerance, maxIter)
	fmt.Printf("Минимум функции достигается в точке (x = %.4f, y = %.4f) со значением f(x,y) = %.4f\n", xMin, yMin, f(xMin, yMin))
	fmt.Println("=================================================================================================\n")

	fmt.Println("Метод Наискорейшего Спуска")
	fmt.Println("--------------------------------------------------------------------------------------------------")
	xMin, yMin = fastestDescent(x0, y0, tolerance, maxIter)
	fmt.Printf("Минимум функции достигается в точке (x = %.4f, y = %.4f) со значением f(x,y) = %.4f\n", xMin, yMin, f(xMin, yMin))
	fmt.Println("=================================================================================================\n")
}
