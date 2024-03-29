#pragma once

template <typename type_data>
class data3d_sycl {
	type_data* data;
	int n, m, k, delta_x, delta_y, delta_z;
public:
	data3d_sycl(sycl::queue q_, int n, int m, int k, int delta_x, int delta_y, int delta_z) {
		size_t array_size = (n + 2 * delta_x + 2) * (m + 2 * delta_y + 2) * (k + 2 * delta_z + 2);
		this->data = malloc_shared<type_data>(array_size, q_);

		this->n = n; this->m = m; this->k = k;
		this->delta_x = delta_x; this->delta_y = delta_y; this->delta_z = delta_z;
	}
	inline type_data& operator()(int i, int j, int l) {
		return data[i * (m + 2 * delta_y + 2) * (k + 2 * delta_z + 2) + j * (k + 2 * delta_z + 2) + l];
	}
	type_data* ReturnData() {
		return data;
	}
};

template <typename type_data>
struct Fields {
	data3d_sycl<type_data> Ex, Ey, Ez, Bx, By, Bz;

	Fields(sycl::queue q, int n, int m, int k, int delta_x, int delta_y, int delta_z) :
		Ex(q, n, m, k, delta_x, delta_y, delta_z),
		Ey(q, n, m, k, delta_x, delta_y, delta_z),
		Ez(q, n, m, k, delta_x, delta_y, delta_z),
		Bx(q, n, m, k, delta_x, delta_y, delta_z),
		By(q, n, m, k, delta_x, delta_y, delta_z),
		Bz(q, n, m, k, delta_x, delta_y, delta_z) {}

	void Delete(sycl::queue q) {
		free(this->Ex.ReturnData(), q);
		free(this->Ey.ReturnData(), q);
		free(this->Ez.ReturnData(), q);
		free(this->Bx.ReturnData(), q);
		free(this->By.ReturnData(), q);
		free(this->Bz.ReturnData(), q);
	}
};

template <typename type_data>
struct SplitFields {
	data3d_sycl<type_data> Exy, Exz, Eyx, Eyz, Ezx, Ezy,
		Bxy, Bxz, Byx, Byz, Bzx, Bzy;
	SplitFields(sycl::queue q, int n, int m, int k, int delta_x, int delta_y, int delta_z) :
		Exy(q, n, m, k, delta_x, delta_y, delta_z),
		Exz(q, n, m, k, delta_x, delta_y, delta_z),
		Eyx(q, n, m, k, delta_x, delta_y, delta_z),
		Eyz(q, n, m, k, delta_x, delta_y, delta_z),
		Ezx(q, n, m, k, delta_x, delta_y, delta_z),
		Ezy(q, n, m, k, delta_x, delta_y, delta_z),

		Bxy(q, n, m, k, delta_x, delta_y, delta_z),
		Bxz(q, n, m, k, delta_x, delta_y, delta_z),
		Byx(q, n, m, k, delta_x, delta_y, delta_z),
		Byz(q, n, m, k, delta_x, delta_y, delta_z),
		Bzx(q, n, m, k, delta_x, delta_y, delta_z),
		Bzy(q, n, m, k, delta_x, delta_y, delta_z) {}

	void Delete(sycl::queue q) {
		free(this->Exy.ReturnData(), q);
		free(this->Exz.ReturnData(), q);
		free(this->Eyx.ReturnData(), q);
		free(this->Eyz.ReturnData(), q);
		free(this->Ezx.ReturnData(), q);
		free(this->Ezy.ReturnData(), q);
		free(this->Bxy.ReturnData(), q);
		free(this->Bxz.ReturnData(), q);
		free(this->Byx.ReturnData(), q);
		free(this->Byz.ReturnData(), q);
		free(this->Bzx.ReturnData(), q);
		free(this->Bzy.ReturnData(), q);
	}
};

template <typename type_data>
struct Sigma {
	data3d_sycl<type_data> sigmaEx, sigmaEy, sigmaEz, sigmaBx, sigmaBy, sigmaBz;
	Sigma(sycl::queue q, int n, int m, int k, int delta_x, int delta_y, int delta_z) :
		sigmaEx(q, n, m, k, delta_x, delta_y, delta_z),
		sigmaEy(q, n, m, k, delta_x, delta_y, delta_z),
		sigmaEz(q, n, m, k, delta_x, delta_y, delta_z),
		sigmaBx(q, n, m, k, delta_x, delta_y, delta_z),
		sigmaBy(q, n, m, k, delta_x, delta_y, delta_z),
		sigmaBz(q, n, m, k, delta_x, delta_y, delta_z) {}

	void Delete(sycl::queue q) {
		free(this->sigmaEx.ReturnData(), q);
		free(this->sigmaEy.ReturnData(), q);
		free(this->sigmaEz.ReturnData(), q);
		free(this->sigmaBx.ReturnData(), q);
		free(this->sigmaBy.ReturnData(), q);
		free(this->sigmaBz.ReturnData(), q);
	}
};
template <typename type_data>
struct Coefficient {
	data3d_sycl<type_data> Exy1, Exz1, Eyx1, Eyz1, Ezx1, Ezy1, Bxy1, Bxz1, Byx1, Byz1, Bzx1, Bzy1,
		Exy2, Exz2, Eyx2, Eyz2, Ezx2, Ezy2, Bxy2, Bxz2, Byx2, Byz2, Bzx2, Bzy2;

	Coefficient(sycl::queue q, int n, int m, int k, int delta_x, int delta_y, int delta_z) :
		Exy1(q, n, m, k, delta_x, delta_y, delta_z),
		Exz1(q, n, m, k, delta_x, delta_y, delta_z),
		Eyx1(q, n, m, k, delta_x, delta_y, delta_z),
		Eyz1(q, n, m, k, delta_x, delta_y, delta_z),
		Ezx1(q, n, m, k, delta_x, delta_y, delta_z),
		Ezy1(q, n, m, k, delta_x, delta_y, delta_z),

		Bxy1(q, n, m, k, delta_x, delta_y, delta_z),
		Bxz1(q, n, m, k, delta_x, delta_y, delta_z),
		Byx1(q, n, m, k, delta_x, delta_y, delta_z),
		Byz1(q, n, m, k, delta_x, delta_y, delta_z),
		Bzx1(q, n, m, k, delta_x, delta_y, delta_z),
		Bzy1(q, n, m, k, delta_x, delta_y, delta_z),

		Exy2(q, n, m, k, delta_x, delta_y, delta_z),
		Exz2(q, n, m, k, delta_x, delta_y, delta_z),
		Eyx2(q, n, m, k, delta_x, delta_y, delta_z),
		Eyz2(q, n, m, k, delta_x, delta_y, delta_z),
		Ezx2(q, n, m, k, delta_x, delta_y, delta_z),
		Ezy2(q, n, m, k, delta_x, delta_y, delta_z),

		Bxy2(q, n, m, k, delta_x, delta_y, delta_z),
		Bxz2(q, n, m, k, delta_x, delta_y, delta_z),
		Byx2(q, n, m, k, delta_x, delta_y, delta_z),
		Byz2(q, n, m, k, delta_x, delta_y, delta_z),
		Bzx2(q, n, m, k, delta_x, delta_y, delta_z),
		Bzy2(q, n, m, k, delta_x, delta_y, delta_z) {}

	void Delete(sycl::queue q) {
		free(this->Exy1.ReturnData(), q);
		free(this->Exz1.ReturnData(), q);
		free(this->Eyx1.ReturnData(), q);
		free(this->Eyz1.ReturnData(), q);
		free(this->Ezx1.ReturnData(), q);
		free(this->Ezy1.ReturnData(), q);
		free(this->Bxy1.ReturnData(), q);
		free(this->Bxz1.ReturnData(), q);
		free(this->Byx1.ReturnData(), q);
		free(this->Byz1.ReturnData(), q);
		free(this->Bzx1.ReturnData(), q);
		free(this->Bzy1.ReturnData(), q);

		free(this->Exy2.ReturnData(), q);
		free(this->Exz2.ReturnData(), q);
		free(this->Eyx2.ReturnData(), q);
		free(this->Eyz2.ReturnData(), q);
		free(this->Ezx2.ReturnData(), q);
		free(this->Ezy2.ReturnData(), q);
		free(this->Bxy2.ReturnData(), q);
		free(this->Bxz2.ReturnData(), q);
		free(this->Byx2.ReturnData(), q);
		free(this->Byz2.ReturnData(), q);
		free(this->Bzx2.ReturnData(), q);
		free(this->Bzy2.ReturnData(), q);
	}
};