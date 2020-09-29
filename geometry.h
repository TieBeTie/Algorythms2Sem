#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>
#include <assert.h>

using std::cin;
using std::cout;
using std::vector;
using std::pow;
using std::max;
using std::min;
using std::pair;
using std::swap;

const double EPS = 1.e-6;
const double PI = 3.1415926535897932384626433832795028841971;

struct Point {
    Point() {}
    Point(double x, double y) {
        this->x = x;
        this->y = y;
    }
    double x, y;
};
bool operator==(const Point& p1, const Point& p2) {
    return abs(p1.x - p2.x) < EPS && abs(p1.y - p2.y) < EPS;
}
bool operator!=(const Point& p1, const Point& p2) {
    return !(p1 == p2);
}
Point operator+ (const Point& p1, const Point& p2) {
    return Point(p1.x + p2.x, p1.y + p2.y);
}
Point operator* (const double k, const Point& p) {
    return Point(k * p.x, k * p.y);
}
Point operator/ (const Point& p, const double k) {
    return Point(p.x / k, p.y / k);
}
bool IsEqual(const double a, const double b) {
    return abs(a - b) < EPS;
}

double Length(const Point& p1, const Point& p2) {
    return sqrt((p1.x - p2.x) * (p1.x - p2.x)
        + (p1.y - p2.y) * (p1.y - p2.y));
}

struct Line {
    Line(double a, double b, double c) {
        this->a = a;
        this->b = b;
        this->c = c;
    }
    Line(Point p1, Point p2) {
        a = p1.y - p2.y;
        b = p2.x - p1.x;
        c = p1.x * p2.y - p1.y * p2.x;
    }
    Line(double k, double s) {
        a = k;
        b = -1;
        c = s;
    }
    Line(Point p, double k) {
        a = k;
        b = -1;
        c = p.y - k * p.x;
    }
    double a, b, c;
};
bool operator==(const Line& l1, const Line& l2) {
    double k;
    if (!IsEqual(l1.a, 0) && IsEqual(l2.a, 0) ||
        IsEqual(l1.a, 0) && !IsEqual(l2.a, 0)) {
        return false;
    }
    else if (!IsEqual(l1.a, 0) && !IsEqual(l2.a, 0)) {
        k = l1.a / l2.a;
    }
    
    if (!IsEqual(l1.b, 0) && IsEqual(l2.b, 0) ||
        IsEqual(l1.b, 0) && !IsEqual(l2.b, 0)) {
        return false;
    }
    else if (!IsEqual(l1.b, 0) && !IsEqual(l2.b, 0)) {
        k = l1.b / l2.b;
    }

    if (!IsEqual(l1.c, 0) && IsEqual(l2.c, 0) ||
        IsEqual(l1.c, 0) && !IsEqual(l2.c, 0)) {
        return false;
    }
    else if (!IsEqual(l1.a, 0) && !IsEqual(l2.a, 0)) {
        k = l1.c / l2.c;
    }


    return IsEqual(l1.a / l2.a, k) &&
        IsEqual(l1.b / l2.b, k) &&
        IsEqual(l1.c / l2.c, k);
}
bool operator!=(const Line& l1, const Line& l2) {
    return !(l1 == l2);
}
Point Intersection(const Line& l1, const Line& l2) {
    return Point(-((l1.c * l2.b - l2.c * l1.b) / (l1.a * l2.b - l2.a * l1.b)),
        -((l1.a * l2.c - l2.a * l1.c) / (l1.a * l2.b - l2.a * l1.b)));
}

struct Vector {
    double x, y;
    Vector(Point A, Point B) {
        x = B.x - A.x;
        y = B.y - A.y;
    }
    Vector(double x, double y) {
        this->x = x;
        this->y = y;
    }
};
Vector operator* (const double k, const Vector& v) {
    return Vector(k * v.x, k * v.y);
}
Vector operator+ (const Vector& v1, const Vector& v2) {
    return Vector(v1.x + v2.x, v1.y + v2.y);
}
double operator* (const Vector& v1, const Vector& v2) {
    return (v1.x * v2.x + v1.y * v2.y) /
        (sqrt(v1.x * v1.x + v1.y * v1.y) * sqrt(v2.x * v2.x + v2.y * v2.y));
}
double operator==(const Vector& v1, const Vector& v2) {
    return   abs(v1.x - v2.x) < EPS && abs(v1.y - v2.y) < EPS;
}
double operator!=(const Vector& v1, const Vector& v2) {
    return !(abs(v1.x - v2.x) < EPS && abs(v1.y - v2.y) < EPS);
}

class Shape {
public:
    virtual ~Shape() = default;
    virtual double perimeter() const = 0;
    virtual double area() const = 0;
    virtual bool isCongruentTo(const Shape& another) const = 0;
    virtual bool isSimilarTo(const Shape& another) const = 0;
    virtual bool containsPoint(Point point) const = 0;
    virtual void rotate(Point center, double angle) = 0;
    virtual void reflex(Point center) = 0;
    virtual void reflex(Line axis) = 0;
    virtual void scale(Point center, double coefficient) = 0;
};

template <class Base, class T>
bool IsInstance(const T* arg) {
    return dynamic_cast <const Base*> (arg);
}

Point RotatePoint(Point& p, Point const center, double angle) {
    angle = angle * PI / 180;
    double x = center.x + (p.x - center.x) * cos(angle) - (p.y - center.y) * sin(angle);
    double y = center.y + (p.x - center.x) * sin(angle) + (p.y - center.y) * cos(angle);
    return Point(x, y);
}


Point ScalePoint(Point& p, Point const center, double coefficient) {
    double x = center.x + (p.x - center.x) * coefficient;
    double y = center.y + (p.y - center.y) * coefficient;
    return Point(x, y);
}


Line Perpendicular(Point const p, Line const l) {
    return Line(-l.b, l.a, -(l.a * p.y) + (l.b * p.x));
}

double Angle(Vector const a, Vector const b) {
    return acos(a * b);
}

class Polygon : public Shape {
public:
    Polygon() {}
    template <typename ... Points>
    Polygon(const Points& ... points) {
        p = vector <Point>{ points... };
        count = p.size();
    }
    Polygon(vector <Point>& points) {
        p = points;
        count = p.size();
    }
    double perimeter() const override {
        double sum = 0;
        for (size_t i = 0; i < count; ++i) {
            sum += Length(p[i], p[(i + 1) % count]);
        }
        return sum;
    }
    double area() const override {
        double sum = 0.0;
        for (size_t i = 0; i < count; ++i) {
            sum += (p[i].x + p[(i + 1) % count].x)
                * (p[i].y - p[(i + 1) % count].y);
        }
        return abs(sum) / 2.0;
    }
    bool operator==(const Shape& another) const {
        if (!IsInstance <Polygon, Shape>(&another)) {
            return false;
        }

        const Polygon* other = dynamic_cast<const Polygon*>(&another);
        if (count != other->count) {
            return false;
        }

        long long n = count;
        long long same = 0;
        for (long long i = 0; i < n; ++i) {
            if (p[0] == other->p[i]) {
                same = i;
            }
        }

        bool equal = true;
        for (long long i = 0; i < n; ++i) {
            if (p[i] != other->p[(same + i) % n]) {
                equal = false;
                break;
            }
        }
        if (equal) {
            return true;
        }

        equal = true;
        for (long long i = 0; i < n; ++i) {
            if (p[i] != other->p[(same - i + n) % n]) {
                equal = false;
                break;
            }
        }

        return equal;
    }
    bool operator!=(const Shape& another) const {
        return !(*this == another);
    }
    bool isSimilarTo(const Shape& another) const override {

        if (!IsInstance <Polygon, Shape>(&another)) {
            return false;
        }

        const Polygon* other = dynamic_cast<const Polygon*>(&another);
        if (count != other->count) {
            return false;
        }


        long long n = count;
        bool congruent = false;
        double k = sqrt(area() / other->area());

        for (long long same = 0; same < n; ++same) {
            if (!IsEqual(
                Length(p[0], p[1])
                /
                Length(other->p[same], other->p[(same + 1) % n]), k)) {
                continue;
            }

            congruent = true;
            for (long long j = 0; j < n; ++j) {

                if (!(IsEqual(

                    Length(p[j], p[(1 + j) % n])
                    /
                    Length(
                        other->p[(same + j) % n],
                        other->p[(same + 1 + j) % n]),

                    k))

                    ||

                    !(IsEqual(
                        Vector(p[(1 + j) % n], p[j])
                        *
                        Vector(p[(1 + j) % n], p[(2 + j) % n])
                        ,
                        Vector(
                            other->p[(same + 1 + j) % n],
                            other->p[(same + j) % n])
                        *
                        Vector(
                            other->p[(same + 1 + j) % n],
                            other->p[(same + 2 + j) % n])))
                    ) {
                    congruent = false;
                    break;
                }
            }

            if (congruent) {
                return true;
            }

            congruent = true;
            for (long long j = 0; j < n; ++j) {

                if (!(IsEqual(

                    Length(p[j], p[(1 + j) % n])
                    /
                    Length(
                        other->p[(same + j) % n],
                        other->p[(same + 1 + j) % n]),

                    k))
                    ||
                    !(IsEqual(
                        Vector(p[(1 + j) % n], p[j]) *
                        Vector(p[(1 + j) % n], p[(2 + j) % n]),

                        Vector(other->p[(same - j + n) % n], other->p[(same + 1 - j + n) % n]) *
                        Vector(other->p[(same - j + n) % n], other->p[(same - 1 - j + 2 * n) % n])))
                    ) {
                    congruent = false;
                    break;
                }
            }
            if (congruent) {
                return true;
            }
        }
        return false;
    }
    bool containsPoint(Point point) const override {
        bool contain = false;
        for (size_t i = 0; i < count; ++i) {
            if (point.x == p[i].x && point.y == p[i].y
                || point.x == p[(i + 1) % count].x && point.y == p[(i + 1) % count].y) {
                return true;
            }
            Line linus(p[i], p[(i + 1) % count]);
            if (IsEqual(linus.a * point.x + linus.b * point.y + linus.c, 0) &&
                ((p[i].y > point.y&& point.y > p[(i + 1) % count].y) ||
                (p[i].y < point.y && point.y < p[(i + 1) % count].y)) &&
                    (p[i].x > point.x&& point.x > p[(i + 1) % count].x) ||
                (p[i].x < point.x && point.x < p[(i + 1) % count].x)) {
                return true;
            }

            if (IsEqual(p[i].y, p[(i + 1) % count].y)) {
                continue;
            }

            if (IsEqual(point.y, min(p[i].y, p[(i + 1) % count].y))) {
                continue;
            }

            if (IsEqual(point.y, max(p[i].y, p[(i + 1) % count].y)) &&
                point.x < min(p[i].x, p[(i + 1) % count].x)) {
                contain = !contain;
                continue;
            }


            if (IsEqual(point.y, max(p[i].y, p[(i + 1) % count].y))) {
                continue;
            }

            Line line(p[i], p[(i + 1) % count]);
            if (point.x < -(line.b * point.y + line.c) / line.a &&
                ((p[i].y > point.y&& point.y > p[(i + 1) % count].y) ||
                (p[i].y < point.y && point.y < p[(i + 1) % count].y))) {
                contain = !contain;
            }
        }
        return contain;
    }
    bool isCongruentTo(const Shape& another) const override {
        if (!IsInstance <Polygon, Shape>(&another)) {
            return false;
        }

        const Polygon* other = dynamic_cast<const Polygon*>(&another);
        if (count != other->count) {
            return false;
        }

        if (!(IsEqual(area(), other->area())) ||
            !(IsEqual(perimeter(), other->perimeter()))) {
            return false;
        }

        long long n = count;
        bool congruent;

        for (long long same = 0; same < n; ++same) {
            if (!IsEqual(Length(p[0], p[1]),
                Length(other->p[same], other->p[(same + 1) % n]))) {
                continue;
            }

            congruent = true;
            for (long long j = 0; j < n; ++j) {

                if (!(IsEqual(
                    Length(p[j], p[(1 + j) % n]),
                    Length(
                        other->p[(same + j) % n],
                        other->p[(same + 1 + j) % n])))
                    ||

                    !(IsEqual(
                        Vector(p[(1 + j) % n], p[j])
                        *
                        Vector(p[(1 + j) % n], p[(2 + j) % n])
                        ,
                        Vector(
                            other->p[(same + 1 + j) % n],
                            other->p[(same + j) % n])
                        *
                        Vector(
                            other->p[(same + 1 + j) % n],
                            other->p[(same + 2 + j) % n])))
                    ) {
                    congruent = false;
                    break;
                }
            }

            if (congruent) {
                return true;
            }

            congruent = true;
            for (long long j = 0; j < n; ++j) {

                if (!(IsEqual(
                    Length(p[j], p[(1 + j) % n]),
                    Length(
                        other->p[(same - j + n) % n],
                        other->p[(same + 1 - j + n) % n])))
                    ||
                    !(IsEqual(
                        Vector(p[(1 + j) % n], p[j]) *
                        Vector(p[(1 + j) % n], p[(2 + j) % n]),

                        Vector(other->p[(same - j + n) % n], other->p[(same + 1 - j + n) % n]) *
                        Vector(other->p[(same - j + n) % n], other->p[(same - 1 - j + 2 * n) % n])))
                    ) {
                    congruent = false;
                    break;
                }
            }
            if (congruent) {
                return true;
            }

        }
        return false;
    }
    void rotate(Point center, double angle) override {
        for (size_t i = 0; i < count; ++i) {
            p[i] = RotatePoint(p[i], center, angle);
        }
    }
    void reflex(Point center) override {

        for (size_t i = 0; i < count; ++i) {
            p[i].x = -p[i].x + center.x;
            p[i].y = -p[i].y + center.y;
        }
    }
    void reflex(Line axis) override {

        Point tmp;
        for (size_t i = 0; i < count; ++i) {
            Line perpen(-axis.b, axis.a, -(axis.a * p[i].y) + (axis.b * p[i].x));
            tmp = Intersection(axis, perpen);
            p[i].x = 2 * tmp.x - p[i].x;
            p[i].y = 2 * tmp.y - p[i].y;
        }

    }
    void scale(Point center, double coefficient) override {
        for (size_t i = 0; i < count; ++i) {
            p[i] = ScalePoint(p[i], center, coefficient);
        }
    }
    bool isConvex() {
        int rotater = -1;
        for (size_t i = 0; i < count; ++i) {
            Line line(p[i], p[(i + 1) % count]);
            Point point = p[(i + 2) % count];
            if (line.a * point.x + line.b * point.y + line.c > 0) {
                rotater = 1;
                break;
            }
        }
        for (size_t i = 0; i < count; ++i) {
            Line line(p[i], p[(i + 1) % count]);
            Point point = p[(i + 2) % count];
            if (rotater * (line.a * point.x + line.b * point.y + line.c) < 0) {
                return false;
            }
        }
        return true;
    }
    size_t verticesCount() {
        return count;
    }
    const vector<Point>& getVertices() {
        return p;
    }
protected:
    vector <Point> p;
    size_t count;
};

class Ellipse : public Shape {
public:
    Ellipse() {}
    Ellipse(Point f1, Point f2, double dist)  {
        this->f1 = f1;
        this->f2 = f2;
        this->dist = dist;
        c = Length(f1, f2) / 2.0;
        a = dist / 2.0;
        e = c / a;
        b = a * sqrt(1 - e * e);
        p = a * (1 - e * e);
    }
    double perimeter() const override {
        return PI * (3 * (a + b) - sqrt((3 * a + b) * (a + 3 * b)));
    }
    double area() const override {
        return PI * a * b;
    }
    bool operator==(const Shape& another) const {
        if (!IsInstance <Ellipse, Shape>(&another)) {
            return false;
        }

        const Ellipse* another_ell = dynamic_cast<const Ellipse*>(&another);
        return IsEqual(dist, another_ell->dist) &&
            (((f1 == another_ell->f1) && (f2 == another_ell->f2)) ||
            ((f1 == another_ell->f2) && (f2 == another_ell->f1)));
    }
    bool operator!=(const Shape& another) const {
        return !(*this == another);
    }
    bool isCongruentTo(const Shape& another) const override {
        if (!IsInstance <Ellipse, Shape>(&another)) {
            return false;
        }
        const Ellipse* another_ell = dynamic_cast<const Ellipse*>(&another);
        return IsEqual(a, another_ell->a) && IsEqual(b, another_ell->b);
    }
    bool isSimilarTo(const Shape& another) const override {
        if (!IsInstance <Ellipse, Shape>(&another)) {
            return false;
        }

        const Ellipse* another_ell = dynamic_cast<const Ellipse*>(&another);
        return IsEqual(a / b, another_ell->a / another_ell->b);
    }
    bool containsPoint(Point point) const override {
        return Length(f1, point) + Length(f2, point) < dist + EPS;
    }
    void rotate(Point center, double angle) override {
        f1 = RotatePoint(f1, center, angle);
        f2 = RotatePoint(f2, center, angle);
    }
    void reflex(Point center) override {
        f1.x = -f1.x + center.x;
        f1.y = -f1.y + center.y;
        f2.x = -f2.x + center.x;
        f2.y = -f2.y + center.y;
    }
    void reflex(Line axis) override {
        Point tmp;
        Line perpen_f1(-axis.b, axis.a, -(axis.a * f1.y) + (axis.b * f1.x));
        tmp = Intersection(axis, perpen_f1);
        f1.x = 2 * tmp.x - f1.x;
        f1.y = 2 * tmp.y - f1.y;

        Line perpen_f2(-axis.b, axis.a, -(axis.a * f2.y) + (axis.b * f2.x));
        tmp = Intersection(axis, perpen_f2);
        f2.x = 2 * tmp.x - f2.x;
        f2.y = 2 * tmp.y - f2.y;
    }
    void scale(Point center, double coefficient) override {
        f1 = ScalePoint(f1, center, coefficient);
        f2 = ScalePoint(f2, center, coefficient);
        dist *= coefficient;
    }
    pair <Point, Point> focuses() const  {
        return pair <Point, Point>(f1, f2);
    }
    pair <Line, Line> directrices() const {
        return pair <Line, Line>(Line(1, 0, -p / (e * (1 + e))),
            Line(1, 0, p / (e * (1 + e))));
    }
    double eccentricity() const {
        return e;
    }
    Point center() const  {
        return Point((f1.x + f2.x) / 2.0, (f1.y + f2.y) / 2.0);
    }
protected:
    Point f1, f2;
    double a, b, c, e, p, dist;
};

class Circle : public Ellipse {
public:
    Circle() = default;
    Circle(Point center, double radius) : Ellipse(center, center, radius * 2) {}
    double radius() const  {
        return dist / 2.0;
    }
    Point center() const {
        return f1;
    }
};

class Rectangle : public Polygon {
public:
    Rectangle() {}
    Rectangle(Point p1, Point p2, double ratio) : 
        Polygon(p1, Point(0, 0), p2, Point(0, 0)) {
        if (ratio > 1) {
            ratio = 1 / ratio;
        }

        double angle = atan(ratio) * 180 / PI;
        double turn = 180 - 2 * angle;

        p[1] = p[0];
        p[0] = RotatePoint(p[0], Point((p[0] + p[2]) / 2.0), turn);
        swap(p[0], p[1]);

        p[3] = p[2];
        p[2] = RotatePoint(p[2], Point((p[0] + p[2]) / 2.0), turn);
        swap(p[2], p[3]);
    }
    Point center() const {
        return Point((p[0] + p[2]) / 2.0);
    }
    pair <Line, Line> diagonals() const {
        return pair <Line, Line>(Line(p[0], p[2]), Line(p[1], p[3]));
    }
};

class Square : public Rectangle {
public:
    Square(Point p1, Point p2) : Rectangle(p1, p2, 1) {}
    Circle circumscribedCircle() const {
        return Circle(center(), Length(center(), p[0]));
    }
    Circle inscribedCircle() const {
        return Circle(center(), Length(p[0], p[1]) / 2.0);
    }
};

class Triangle : public Polygon {
public:
    Triangle() = default;
    template <typename ... Points>
    Triangle(const Points& ... points) : Polygon(points...) {
        p = vector <Point>{ points... };
        count = p.size();
    }
    Circle circumscribedCircle() const {
        Point middle_a = (p[0] + p[1]) / 2.0;
        Line line_a(p[0], p[1]);
        Line perpen_a = Perpendicular(middle_a, line_a);

        Point middle_b = (p[1] + p[2]) / 2.0;
        Line line_b(p[1], p[2]);
        Line perpen_b = Perpendicular(middle_b, line_b);
        
        return Circle(Intersection(perpen_a, perpen_b), 
            (Length(p[0], p[1]) * Length(p[1], p[2]) * Length(p[2], p[0])) 
            / (4.0 * area()));
    }
    Circle inscribedCircle() const {

        double A = Length(p[1], p[2]);
        double B = Length(p[2], p[0]);
        double C = Length(p[0], p[1]);

        Point incenter  = (A * p[0] + B * p[1] + C * p[2])
            / (A + B + C);

        double per = perimeter() / 2.0;
        double r = area() / per;
        return Circle(incenter, r);
    }
    Point centroid() const {
        return (p[0] + p[1] + p[2]) / 3.0;
    }
    Point orthocenter() const {
       Line perpen_a = Perpendicular(p[0], Line(p[1], p[2]));
       Line perpen_b = Perpendicular(p[1], Line(p[2], p[0]));
       return Intersection(perpen_a, perpen_b);
    }
    Line EulerLine() const {
        return Line(circumscribedCircle().center(), orthocenter());
    }
    Circle ninePointsCircle() const {
        Point nine_center = (circumscribedCircle().center() + orthocenter()) / 2.0;
        return Circle(nine_center, circumscribedCircle().radius() / 2.0);
    }
};
