#include "Render.h"

#include <sstream>
#include "Render.h"

#include <sstream>
#include <iostream>
#include <stdlib.h>
#include <windows.h>
#include <GL\GL.h>
#include <GL\GLU.h>

#include "MyOGL.h"

#include "Camera.h"
#include "Light.h"
#include "Primitives.h"

#include "GUItextRectangle.h"
#include <random>
bool textureMode = true;
bool lightMode = true;

//класс для настройки камеры
class CustomCamera : public Camera
{
public:
	//дистанция камеры
	double camDist;
	//углы поворота камеры
	double fi1, fi2;

	
	//значния масеры по умолчанию
	CustomCamera()
	{
		camDist = 15;
		fi1 = 1;
		fi2 = 1;
	}

	
	//считает позицию камеры, исходя из углов поворота, вызывается движком
	void SetUpCamera()
	{
		//отвечает за поворот камеры мышкой
		lookPoint.setCoords(0, 0, 0);

		pos.setCoords(camDist*cos(fi2)*cos(fi1),
			camDist*cos(fi2)*sin(fi1),
			camDist*sin(fi2));

		if (cos(fi2) <= 0)
			normal.setCoords(0, 0, -1);
		else
			normal.setCoords(0, 0, 1);

		LookAt();
	}

	void CustomCamera::LookAt()
	{
		//функция настройки камеры
		gluLookAt(pos.X(), pos.Y(), pos.Z(), lookPoint.X(), lookPoint.Y(), lookPoint.Z(), normal.X(), normal.Y(), normal.Z());
	}



}  camera;   //создаем объект камеры


//Класс для настройки света
class CustomLight : public Light
{
public:
	CustomLight()
	{
		//начальная позиция света
		pos = Vector3(1, 1, 3);
	}

	
	//рисует сферу и линии под источником света, вызывается движком
	void  DrawLightGhismo()
	{
		glDisable(GL_LIGHTING);

		
		glColor3d(0.9, 0.8, 0);
		Sphere s;
		s.pos = pos;
		s.scale = s.scale*0.08;
		s.Show();
		
		if (OpenGL::isKeyPressed('G'))
		{
			glColor3d(0, 0, 0);
			//линия от источника света до окружности
			glBegin(GL_LINES);
			glVertex3d(pos.X(), pos.Y(), pos.Z());
			glVertex3d(pos.X(), pos.Y(), 0);
			glEnd();

			//рисуем окруность
			Circle c;
			c.pos.setCoords(pos.X(), pos.Y(), 0);
			c.scale = c.scale*1.5;
			c.Show();
		}

	}

	void SetUpLight()
	{
		GLfloat amb[] = { 0.2, 0.2, 0.2, 0 };
		GLfloat dif[] = { 1.0, 1.0, 1.0, 0 };
		GLfloat spec[] = { .7, .7, .7, 0 };
		GLfloat position[] = { pos.X(), pos.Y(), pos.Z(), 1. };

		// параметры источника света
		glLightfv(GL_LIGHT0, GL_POSITION, position);
		// характеристики излучаемого света
		// фоновое освещение (рассеянный свет)
		glLightfv(GL_LIGHT0, GL_AMBIENT, amb);
		// диффузная составляющая света
		glLightfv(GL_LIGHT0, GL_DIFFUSE, dif);
		// зеркально отражаемая составляющая света
		glLightfv(GL_LIGHT0, GL_SPECULAR, spec);

		glEnable(GL_LIGHT0);
	}


} light;  //создаем источник света




//старые координаты мыши
int mouseX = 0, mouseY = 0;

void mouseEvent(OpenGL *ogl, int mX, int mY)
{
	int dx = mouseX - mX;
	int dy = mouseY - mY;
	mouseX = mX;
	mouseY = mY;

	//меняем углы камеры при нажатой левой кнопке мыши
	if (OpenGL::isKeyPressed(VK_RBUTTON))
	{
		camera.fi1 += 0.01*dx;
		camera.fi2 += -0.01*dy;
	}

	
	//двигаем свет по плоскости, в точку где мышь
	if (OpenGL::isKeyPressed('G') && !OpenGL::isKeyPressed(VK_LBUTTON))
	{
		LPPOINT POINT = new tagPOINT();
		GetCursorPos(POINT);
		ScreenToClient(ogl->getHwnd(), POINT);
		POINT->y = ogl->getHeight() - POINT->y;

		Ray r = camera.getLookRay(POINT->x, POINT->y);

		double z = light.pos.Z();

		double k = 0, x = 0, y = 0;
		if (r.direction.Z() == 0)
			k = 0;
		else
			k = (z - r.origin.Z()) / r.direction.Z();

		x = k*r.direction.X() + r.origin.X();
		y = k*r.direction.Y() + r.origin.Y();

		light.pos = Vector3(x, y, z);
	}

	if (OpenGL::isKeyPressed('G') && OpenGL::isKeyPressed(VK_LBUTTON))
	{
		light.pos = light.pos + Vector3(0, 0, 0.02*dy);
	}

	
}

void mouseWheelEvent(OpenGL *ogl, int delta)
{

	if (delta < 0 && camera.camDist <= 1)
		return;
	if (delta > 0 && camera.camDist >= 100)
		return;

	camera.camDist += 0.01*delta;

}

void keyDownEvent(OpenGL *ogl, int key)
{
	if (key == 'L')
	{
		lightMode = !lightMode;
	}

	if (key == 'T')
	{
		textureMode = !textureMode;
	}

	if (key == 'R')
	{
		camera.fi1 = 1;
		camera.fi2 = 1;
		camera.camDist = 15;

		light.pos = Vector3(1, 1, 3);
	}

	if (key == 'F')
	{
		light.pos = camera.pos;
	}
}

void keyUpEvent(OpenGL *ogl, int key)
{
	
}

double fRand(double fMin, double fMax)
{
	double f = (double)rand() / RAND_MAX;
	return fMin + f * (fMax - fMin);
}

GLuint texId;
Vector3 points[4][4];


//выполняется перед первым рендером
void initRender(OpenGL *ogl)
{
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
		{
			points[i][j] = Vector3(i * 5, -(j + 0.3) * 3, fRand(-1.0, 5.0));
		}
	//настройка текстур

	//4 байта на хранение пикселя
	glPixelStorei(GL_UNPACK_ALIGNMENT, 4);

	//настройка режима наложения текстур
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

	//включаем текстуры
	glEnable(GL_TEXTURE_2D);
	

	//массив трехбайтных элементов  (R G B)
	RGBTRIPLE *texarray;

	//массив символов, (высота*ширина*4      4, потомучто   выше, мы указали использовать по 4 байта на пиксель текстуры - R G B A)
	char *texCharArray;
	int texW, texH;
	OpenGL::LoadBMP("texture.bmp", &texW, &texH, &texarray);
	OpenGL::RGBtoChar(texarray, texW, texH, &texCharArray);

	
	
	//генерируем ИД для текстуры
	glGenTextures(1, &texId);
	//биндим айдишник, все что будет происходить с текстурой, будте происходить по этому ИД
	glBindTexture(GL_TEXTURE_2D, texId);

	//загружаем текстуру в видеопямять, в оперативке нам больше  она не нужна
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, texW, texH, 0, GL_RGBA, GL_UNSIGNED_BYTE, texCharArray);

	//отчистка памяти
	free(texCharArray);
	free(texarray);

	//наводим шмон
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST); 
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);


	//камеру и свет привязываем к "движку"
	ogl->mainCamera = &camera;
	ogl->mainLight = &light;

	// нормализация нормалей : их длины будет равна 1
	glEnable(GL_NORMALIZE);

	// устранение ступенчатости для линий
	glEnable(GL_LINE_SMOOTH); 


	//   задать параметры освещения
	//  параметр GL_LIGHT_MODEL_TWO_SIDE - 
	//                0 -  лицевые и изнаночные рисуются одинаково(по умолчанию), 
	//                1 - лицевые и изнаночные обрабатываются разными режимами       
	//                соответственно лицевым и изнаночным свойствам материалов.    
	//  параметр GL_LIGHT_MODEL_AMBIENT - задать фоновое освещение, 
	//                не зависящее от сточников
	// по умолчанию (0.2, 0.2, 0.2, 1.0)

	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 0);

	camera.fi1 = -1.3;
	camera.fi2 = 0.8;
}


void Draw_Cube(double P1[3])
{
	glTranslated(-0.5, 0, 0);

	glBegin(GL_QUADS);

	glColor3d(0.7, 0.7, 0.7);
	glVertex3d(P1[0] - 0.5, P1[1] - 0.5, P1[2] - 0.5);
	glVertex3d(P1[0] + 0.5, P1[1] - 0.5, P1[2] - 0.5);
	glVertex3d(P1[0] + 0.5, P1[1] + 0.5, P1[2] - 0.5);
	glVertex3d(P1[0] - 0.5, P1[1] + 0.5, P1[2] - 0.5);

	glColor3d(0.2, 0.2, 0.2);
	glVertex3d(P1[0] - 0.5, P1[1] - 0.5, P1[2] + 0.5);
	glVertex3d(P1[0] + 0.5, P1[1] - 0.5, P1[2] + 0.5);
	glVertex3d(P1[0] + 0.5, P1[1] + 0.5, P1[2] + 0.5);
	glVertex3d(P1[0] - 0.5, P1[1] + 0.5, P1[2] + 0.5);

	glColor3d(0.3, 0.3, 0.3);
	glVertex3d(P1[0] - 0.5, P1[1] - 0.5, P1[2] - 0.5);
	glVertex3d(P1[0] - 0.5, P1[1] + 0.5, P1[2] - 0.5);
	glVertex3d(P1[0] - 0.5, P1[1] + 0.5, P1[2] + 0.5);
	glVertex3d(P1[0] - 0.5, P1[1] - 0.5, P1[2] + 0.5);

	glColor3d(0.4, 0.4, 0.4);
	glVertex3d(P1[0] + 0.5, P1[1] - 0.5, P1[2] - 0.5);
	glVertex3d(P1[0] + 0.5, P1[1] + 0.5, P1[2] - 0.5);
	glVertex3d(P1[0] + 0.5, P1[1] + 0.5, P1[2] + 0.5);
	glVertex3d(P1[0] + 0.5, P1[1] - 0.5, P1[2] + 0.5);

	glColor3d(0.5, 0.5, 0.5);
	glVertex3d(P1[0] - 0.5, P1[1] - 0.5, P1[2] - 0.5);
	glVertex3d(P1[0] + 0.5, P1[1] - 0.5, P1[2] - 0.5);
	glVertex3d(P1[0] + 0.5, P1[1] - 0.5, P1[2] + 0.5);
	glVertex3d(P1[0] - 0.5, P1[1] - 0.5, P1[2] + 0.5);

	glColor3d(0.6, 0.6, 0.6);
	glVertex3d(P1[0] - 0.5, P1[1] + 0.5, P1[2] - 0.5);
	glVertex3d(P1[0] + 0.5, P1[1] + 0.5, P1[2] - 0.5);
	glVertex3d(P1[0] + 0.5, P1[1] + 0.5, P1[2] + 0.5);
	glVertex3d(P1[0] - 0.5, P1[1] + 0.5, P1[2] + 0.5);

	glColor3d(0.2, 0.2, 0.2);
	glVertex3d(P1[0] - 0.5, P1[1] - 0.5, P1[2] - 0.1);
	glVertex3d(P1[0] - 0.5, P1[1] - 0.5, P1[2] + 0.1);
	glVertex3d(P1[0] - 0.5, P1[1] - 1, P1[2] + 0.1);
	glVertex3d(P1[0] - 0.5, P1[1] - 1, P1[2] - 0.1);

	glColor3d(0.3, 0.3, 0.3);
	glVertex3d(P1[0] + 0.5, P1[1] - 0.5, P1[2] - 0.1);
	glVertex3d(P1[0] + 0.5, P1[1] - 0.5, P1[2] + 0.1);
	glVertex3d(P1[0] + 0.5, P1[1] - 1, P1[2] + 0.1);
	glVertex3d(P1[0] + 0.5, P1[1] - 1, P1[2] - 0.1);

	glColor3d(0.4, 0.4, 0.4);
	glVertex3d(P1[0] - 0.5, P1[1] - 1, P1[2] - 0.1);
	glVertex3d(P1[0] - 0.5, P1[1] - 1, P1[2] + 0.1);
	glVertex3d(P1[0] + 0.5, P1[1] - 1, P1[2] + 0.1);
	glVertex3d(P1[0] + 0.5, P1[1] - 1, P1[2] - 0.1);

	glColor3d(0.5, 0.5, 0.5);
	glVertex3d(P1[0] - 0.5, P1[1] - 0.5, P1[2] - 0.1);
	glVertex3d(P1[0] + 0.5, P1[1] - 0.5, P1[2] - 0.1);
	glVertex3d(P1[0] + 0.5, P1[1] - 1, P1[2] - 0.1);
	glVertex3d(P1[0] - 0.5, P1[1] - 1, P1[2] - 0.1);

	glColor3d(0.6, 0.6, 0.6);
	glVertex3d(P1[0] - 0.5, P1[1] - 0.5, P1[2] + 0.1);
	glVertex3d(P1[0] + 0.5, P1[1] - 0.5, P1[2] + 0.1);
	glVertex3d(P1[0] + 0.5, P1[1] - 1, P1[2] + 0.1);
	glVertex3d(P1[0] - 0.5, P1[1] - 1, P1[2] + 0.1);

	glColor3d(0.2, 0.2, 0.2);
	glVertex3d(P1[0] - 0.5, P1[1] + 0.5, P1[2] - 0.1);
	glVertex3d(P1[0] - 0.5, P1[1] + 0.5, P1[2] + 0.1);
	glVertex3d(P1[0] - 0.5, P1[1] + 1, P1[2] + 0.1);
	glVertex3d(P1[0] - 0.5, P1[1] + 1, P1[2] - 0.1);

	glColor3d(0.3, 0.3, 0.3);
	glVertex3d(P1[0] + 0.5, P1[1] + 0.5, P1[2] - 0.1);
	glVertex3d(P1[0] + 0.5, P1[1] + 0.5, P1[2] + 0.1);
	glVertex3d(P1[0] + 0.5, P1[1] + 1, P1[2] + 0.1);
	glVertex3d(P1[0] + 0.5, P1[1] + 1, P1[2] - 0.1);

	glColor3d(0.4, 0.4, 0.4);
	glVertex3d(P1[0] - 0.5, P1[1] + 1, P1[2] - 0.1);
	glVertex3d(P1[0] - 0.5, P1[1] + 1, P1[2] + 0.1);
	glVertex3d(P1[0] + 0.5, P1[1] + 1, P1[2] + 0.1);
	glVertex3d(P1[0] + 0.5, P1[1] + 1, P1[2] - 0.1);

	glColor3d(0.5, 0.5, 0.5);
	glVertex3d(P1[0] - 0.5, P1[1] + 0.5, P1[2] - 0.1);
	glVertex3d(P1[0] + 0.5, P1[1] + 0.5, P1[2] - 0.1);
	glVertex3d(P1[0] + 0.5, P1[1] + 1, P1[2] - 0.1);
	glVertex3d(P1[0] - 0.5, P1[1] + 1, P1[2] - 0.1);

	glColor3d(0.6, 0.6, 0.6);
	glVertex3d(P1[0] - 0.5, P1[1] + 0.5, P1[2] + 0.1);
	glVertex3d(P1[0] + 0.5, P1[1] + 0.5, P1[2] + 0.1);
	glVertex3d(P1[0] + 0.5, P1[1] + 1, P1[2] + 0.1);
	glVertex3d(P1[0] - 0.5, P1[1] + 1, P1[2] + 0.1);


	glEnd();
}

void Draw_Cube()
{
	double P[] = { 0, 0 ,0 };
	Draw_Cube(P);
}

double f(double p0, double p1, double p2, double p3, double t)
{
	return (1 - t) * (1 - t) * (1 - t) * p0 + 3 * t * (1 - t) * (1 - t) * p1 + 3 * t * t * (1 - t) * p2 + t * t * t * p3; 
}

void Bese(double P1[3], double P2[3], double P3[3], double P4[3])
{
	glBegin(GL_LINES);
	glVertex3dv(P1);
	glVertex3dv(P2);
	glVertex3dv(P2);
	glVertex3dv(P3);
	glVertex3dv(P3);
	glVertex3dv(P4);
	glEnd();

	glLineWidth(3); 
	glBegin(GL_LINE_STRIP);
	for (double t = 0; t <= 1.0001; t += 0.01)
	{
		double P[3];
		P[0] = f(P1[0], P2[0], P3[0], P4[0], t);
		P[1] = f(P1[1], P2[1], P3[1], P4[1], t);
		P[2] = f(P1[2], P2[2], P3[2], P4[2], t);
		glVertex3dv(P); 
	}
	glEnd();
	glLineWidth(1); 
}

Vector3 VecBese(Vector3* points, double t)
{	 
	return Vector3( f(points[0].X(), points[1].X(), points[2].X(), points[3].X(), t),
					f(points[0].Y(), points[1].Y(), points[2].Y(), points[3].Y(), t),
					f(points[0].Z(), points[1].Z(), points[2].Z(), points[3].Z(), t));
}

double f_E(double p1, double p4, double vec, double vec2, double t)
{
	return p1 * (2 * pow(t, 3) - 3 * pow(t, 2) + 1) + p4 * (-2 * pow(t, 3) + 3 * pow(t, 2)) + vec * (pow(t, 3) - 2 * pow(t, 2) + t) + vec2 * (pow(t, 3) - pow(t, 2));
}

void Ermit(double p0[3], double p1[3], double p2[3], double p3[3])
{
	//вектор 1 начинается в p0 и заканчивается p1, вектор 2 начинается в p3 и заканчивается p2

	double vector1[] = { p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2] };

	double vector2[] = { p2[0] - p3[0], p2[1] - p3[1], p2[2] - p3[2] };

	glBegin(GL_LINES);
	glVertex3dv(p0);
	glVertex3dv(p1);

	glVertex3dv(p3);
	glVertex3dv(p2);
	glEnd();

	glPointSize(10);

	glBegin(GL_POINTS);

	glVertex3dv(p1);
	glVertex3dv(p2);
	glEnd();


	glLineWidth(3);
	glBegin(GL_LINE_STRIP);
	for (double t = 0; t <= 1.0001; t += 0.01)
	{
		double P[3];
		P[0] = f_E(p0[0], p3[0], vector1[0], vector2[0], t);
		P[1] = f_E(p0[1], p3[1], vector1[1], vector2[1], t);
		P[2] = f_E(p0[2], p3[2], vector1[2], vector2[2], t);
		glVertex3dv(P);
	}
	glEnd();
	glLineWidth(1); 

}

void BuildErmit()
{
	double P1[] = { -10,2,0 };
	double P2[] = { -9, 8, 7 };
	double P3[] = { 0, 2, 0 };
	double P4[] = { -2, 0, 5 };

	Ermit(P1, P2, P3, P4);
}

double fact(int n)
{
	int r = 1;
	for (int i = 1; i <= n; i++)
	{
		r *= i;
	}
	return r;
}

double Bern(int i, int n, double u)
{
	return 1.0 * fact(n) * pow(u, i) * pow(1 - u, n - i) / fact(i) / fact(n - i);
}

Vector3 Bezier(Vector3* points, int n, int m, double u, double v)
{
	Vector3 res(0, 0, 0);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			double Bi = Bern(i, n - 1, u);
			double Bj = Bern(j, m - 1, v);
			res = res + points[i * m + j] * Bi * Bj;
		}
	}
	return res;
}

void BeseAnim()
{
	static bool flag = false;
	double h = 0.01;
	static double t = h;

	double P1[] = { 0, 0, 0 };
	double P2[] = { 30, 22, 8 };
	double P3[] = { -5, 40, 4 };
	double P4[] = { 15, 0, 2 };

	Bese(P1, P2, P3, P4);

	Vector3 Points[4];

	Points[0] = Vector3(P1[0], P1[1], P1[2]);
	Points[1] = Vector3(P2[0], P2[1], P2[2]);
	Points[2] = Vector3(P3[0], P3[1], P3[2]);
	Points[3] = Vector3(P4[0], P4[1], P4[2]);

	Vector3 CurPos = VecBese(Points, t);
	Vector3 PrevPos = VecBese(Points, t - h);
	Vector3 Dir = (CurPos - PrevPos).normolize();
	Vector3 OrigDir(1, 0, 0);
	Vector3 RotX(Dir.X(), Dir.Y(), 0);
	RotX = RotX.normolize();

	double cosU = OrigDir.scalarProisvedenie(RotX);
	Vector3 VecPr = OrigDir.vectProisvedenie(RotX);

	double sinSign = VecPr.Z() / abs(VecPr.Z());
	double U = acos(cosU) * 180 / M_PI * sinSign;

	double cosZU = Vector3(0, 0, 0).scalarProisvedenie(Dir);
	double ZU = acos(Dir.Z()) * 180.0 / M_PI - 90;

	glPushMatrix();
	glTranslated(CurPos.X(), CurPos.Y(), CurPos.Z());
	glRotated(ZU, 0, 1, 0);
	glRotated(U, 0, 0, 1);
	Draw_Cube();
	glPopMatrix();

	if (t >= 1)
		flag = true;
	if (t <= 0)
		flag = false;

	if (flag == false)
		t += h;
	if (flag == true)
		t -= h;
}

double* Normal(double A[3], double B[3], double C[3], bool fl)
{
	double vector1[] = { B[0] - A[0], B[1] - A[1] ,B[2] - A[2] };

	double vector2[] = { C[0] - A[0], C[1] - A[1],C[2] - A[2] };

	double vector_normali[] = { vector1[1] * vector2[2] - vector2[1] * vector1[2], -vector1[0] * vector2[2] + vector2[0] * vector1[2], vector1[0] * vector2[1] - vector2[0] * vector1[1] };

	double lenght = -sqrt(vector_normali[0] * vector_normali[0] + vector_normali[1] * vector_normali[1] + vector_normali[2] * vector_normali[2]);

	if (fl)
		lenght *= -1;

	vector_normali[0] /= lenght;
	vector_normali[1] /= lenght;
	vector_normali[2] /= lenght;

	return vector_normali;
}

void BeseSurface()
{
	double h = 0.1;

	glColor3d(0.4, 0.4, 0.4);

	const double n = 4;
	const double m = 4;

	glPointSize(10);

	glBegin(GL_POINTS);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++) {
			glVertex3dv(points[i][j].toArray());
		}
	glEnd();

	glBegin(GL_LINES);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++) {
			glVertex3dv(points[i][j].toArray());
			j == m - 1 ? glVertex3dv(points[i][j].toArray()) : glVertex3dv(points[i][j + 1].toArray());
			glVertex3dv(points[i][j].toArray());
			i == n - 1 ? glVertex3dv(points[i][j].toArray()) : glVertex3dv(points[i + 1][j].toArray());
		}
	glEnd();

	glColor3d(0.7, 0.7, 0.7);

	Vector3 old1;
	Vector3 old2;
	Vector3 a, b, c, d;

	glBegin(GL_TRIANGLES);
	for (double u = h; u <= 1; u += h)
	{
		old1 = Bezier((Vector3*)points, n, m, u, 0);
		old2 = Bezier((Vector3*)points, n, m, u - h, 0);
		for (double v = h; v <= 1; v += h)
		{
			a = Bezier((Vector3*)points, n, m, u, v);
			b = Bezier((Vector3*)points, n, m, u - h, v);
			c = old2;
			d = old1;


			glNormal3dv(Normal(a.toArrayNC(), b.toArrayNC(), c.toArrayNC(), false));
			glTexCoord2d(u, v);
			glVertex3dv(a.toArray());
			glTexCoord2d(u - h, v);
			glVertex3dv(b.toArray());
			glTexCoord2d(u - h, v - h);
			glVertex3dv(c.toArray());


			glNormal3dv(Normal(a.toArrayNC(), d.toArrayNC(), c.toArrayNC(), true));
			glTexCoord2d(u, v);
			glVertex3dv(a.toArray());
			glTexCoord2d(u, v - h);
			glVertex3dv(d.toArray());
			glTexCoord2d(u - h, v - h);
			glVertex3dv(c.toArray());

			old1 = a;
			old2 = b;
		}
	}
	glEnd();
}

void Render(OpenGL *ogl)
{
	glDisable(GL_TEXTURE_2D);
	glDisable(GL_LIGHTING);

	glEnable(GL_DEPTH_TEST);
	if (textureMode)
		glEnable(GL_TEXTURE_2D);

	if (lightMode)
		glEnable(GL_LIGHTING);


	//альфаналожение
	//glEnable(GL_BLEND);
	//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


	//настройка материала
	GLfloat amb[] = { 0.2, 0.2, 0.1, 1. };
	GLfloat dif[] = { 0.4, 0.65, 0.5, 1. };
	GLfloat spec[] = { 0.9, 0.8, 0.3, 1. };
	GLfloat sh = 0.1f * 256;


	//фоновая
	glMaterialfv(GL_FRONT, GL_AMBIENT, amb);
	//дифузная
	glMaterialfv(GL_FRONT, GL_DIFFUSE, dif);
	//зеркальная
	glMaterialfv(GL_FRONT, GL_SPECULAR, spec); \
		//размер блика
		glMaterialf(GL_FRONT, GL_SHININESS, sh);

	//чтоб было красиво, без квадратиков (сглаживание освещения)
	glShadeModel(GL_SMOOTH);
	//===================================
	//Прогать тут  

	glBindTexture(GL_TEXTURE_2D, texId);

	BuildErmit();

	BeseAnim();

	BeseSurface();

   //Сообщение вверху экрана

	
	glMatrixMode(GL_PROJECTION);	//Делаем активной матрицу проекций. 
	                                //(всек матричные операции, будут ее видоизменять.)
	glPushMatrix();   //сохраняем текущую матрицу проецирования (которая описывает перспективную проекцию) в стек 				    
	glLoadIdentity();	  //Загружаем единичную матрицу
	glOrtho(0, ogl->getWidth(), 0, ogl->getHeight(), 0, 1);	 //врубаем режим ортогональной проекции

	glMatrixMode(GL_MODELVIEW);		//переключаемся на модел-вью матрицу
	glPushMatrix();			  //сохраняем текущую матрицу в стек (положение камеры, фактически)
	glLoadIdentity();		  //сбрасываем ее в дефолт

	glDisable(GL_LIGHTING);



	GuiTextRectangle rec;		   //классик моего авторства для удобной работы с рендером текста.
	rec.setSize(300, 150);
	rec.setPosition(10, ogl->getHeight() - 150 - 10);


	std::stringstream ss;
	ss << "T - вкл/выкл текстур" << std::endl;
	ss << "L - вкл/выкл освещение" << std::endl;
	ss << "F - Свет из камеры" << std::endl;
	ss << "G - двигать свет по горизонтали" << std::endl;
	ss << "G+ЛКМ двигать свет по вертекали" << std::endl;
	ss << "Коорд. света: (" << light.pos.X() << ", " << light.pos.Y() << ", " << light.pos.Z() << ")" << std::endl;
	ss << "Коорд. камеры: (" << camera.pos.X() << ", " << camera.pos.Y() << ", " << camera.pos.Z() << ")" << std::endl;
	ss << "Параметры камеры: R="  << camera.camDist << ", fi1=" << camera.fi1 << ", fi2=" << camera.fi2 << std::endl;
	
	rec.setText(ss.str().c_str());
	rec.Draw();

	glMatrixMode(GL_PROJECTION);	  //восстанавливаем матрицы проекции и модел-вью обратьно из стека.
	glPopMatrix();


	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
}