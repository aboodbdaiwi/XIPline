from tqdm.auto import tqdm
from functions import optpoly
import numpy as np
import sigpy as sp
import time
import matplotlib.pyplot as plt
plt.style.use('dark_background')


def save_helper(save, x, itr, lst_time, lst_err, idx, obj=None):
    """save_helper

    Args:
        save (String): Location.
        x (Float): Variable to be saved.
        itr (Int): Iteration number.
        lst_time (List): List of time values.
        lst_err (List): List of error values.
        idx (Int): Index.
        obj (Float, optional): Objective. Defaults to None.
    """
    if save == None:
        return

    tp = sp.get_array_module(x)
    np.save("%s/time.npy" % save, np.cumsum(lst_time))
    if idx is None:
        tp.save("%s/iter_%03d.npy" % (save, itr), x)
    else:
        tp.save("%s/iter_%03d.npy" % (save, itr), x[idx])
    if obj is not None:
        np.save("%s/obj.npy" % save, obj)
    if lst_err is not None:
        np.save("%s/err.npy" % save, lst_err)


def calc_perc_err(ref, x, ord=2, auto_normalize=True):
    """calc_perc_err

    Args:
        ref (Array): Reference image.
        x (Array): Variable to be compared.
        ord (int, optional): Norm order. Defaults to 2.
        auto_normalize (bool, optional): Normalize. Defaults to True.

    Returns:
        Float: Percent error.
    """
    dev = sp.get_device(x)
    xp = dev.xp
    with dev:
        if auto_normalize:
            p = ref/xp.linalg.norm(ref.ravel(), ord=ord)
            q = x/xp.linalg.norm(x.ravel(), ord=ord)
            err = xp.linalg.norm((p - q).ravel(), ord=ord)
        else:
            err = xp.linalg.norm((ref - x).ravel(), ord=ord)
            err /= xp.linalg.norm(ref.ravel(), ord=ord)
        err = sp.to_device(err, sp.cpu_device)
    if np.isnan(err) or np.isinf(err):
        return 100
    return 100 * err


def cg(num_iters, ptol, A, b, P=None, lamda=None, ref=None,
       save=None, verbose=True, idx=None, draw_output=False):
    r"""Conjugate Gradient.

    Solves for the following optimization problem.

    .. math::
      \min_x \frac{1}{2} \| A x - b \|_2^2

    Inputs:
      num_iters (Int): Maximum number of iterations.
      ptol (Float): l1-percentage tolerance between iterates.
      A (Linop): Forward model.
      b (Array): Measurement.
      P (None or Linop): Preconditioner.
      ref (None or Array): Reference to compare against.
      lamda (None or Float): l2-regularization.
      save (None or String): If specified, path to save iterations and
                             timings.
      verbose (Bool): Print information.
      idx (Tuple): If passed, slice iterates before saving.
      draw_output (Bool): Plot x after each iteration.

    Returns:
      x (Array): Reconstruction.
    """
    device = sp.get_device(b)
    xp = device.xp

    def calculate_objective(x):
        obj = sp.to_device(
            0.5 * xp.linalg.norm(A * x - b)**2, sp.cpu_device)
        if lamda:
            obj += sp.to_device(lamda * xp.linalg.norm(x)
                                ** 2, sp.cpu_device)
        return obj

    with device:

        lst_time = []
        lst_err = None
        objective = []
        tolerance = []
        if type(ref) != type(None):
            ref = sp.to_device(ref, device)
            lst_err = ([], [])

        AHA = A.N
        if lamda is not None:
            AHA = AHA + lamda * sp.linop.Identity(A.ishape)

        if P == None:
            P = sp.linop.Identity(A.ishape)

        # Set-up time.
        start_time = time.perf_counter()
        AHb = A.H(b)
        x = xp.zeros_like(AHb)
        r = AHb - AHA(x)
        z = P(r.copy())
        p = z.copy()
        rzold = xp.real(xp.vdot(r, z))
        end_time = time.perf_counter()

        lst_time.append(end_time - start_time)
        if lst_err != None:
            lst_err[0].append(calc_perc_err(ref, x, ord=1))
            lst_err[1].append(calc_perc_err(ref, x, ord=2))

        save_helper(save, x, 0, lst_time, lst_err, idx)

        if verbose:
            pbar = tqdm(total=num_iters, desc="CG",
                        leave=True)

        if draw_output:
            image_fig = plt.figure(figsize=(6, 13.5), dpi=100)

        objective.append(calculate_objective(x))
        for k in range(0, num_iters):
            x_old = x.copy()
            start_time = time.perf_counter()

            Ap = AHA(p)
            pAp = xp.real(xp.vdot(p, Ap)).item()
            if pAp > 0:
                alpha = rzold / pAp
                sp.axpy(x, alpha, p)
                sp.axpy(r, -alpha, Ap)
                z = P(r.copy())
                rznew = xp.real(xp.vdot(r, z))
                beta = rznew / rzold
                sp.xpay(p, beta, z)
                rzold = rznew

                if draw_output:
                    try:
                        image_slice = x[:, :, np.shape(x)[2]//2]  # 3D
                    except:
                        image_slice = x  # 2D
                    image_slice_cpu = np.rot90(
                        abs(sp.to_device(image_slice, -1)), 3)

                    if num_iters > 144:
                        nrows = 18
                        ncols = 8
                    else:
                        nrows = int(np.floor(np.sqrt(num_iters)))
                        ncols = int(np.ceil(np.sqrt(num_iters)))
                    if k < 144:
                        ax = image_fig.add_subplot(nrows, ncols, k + 1)
                        ax.axis('off')
                        ax.imshow(image_slice_cpu, cmap="gray")
                        ax.set_title('Iter %d' %
                                     (k + 1), fontsize='x-small', y=0.92)

            objective.append(calculate_objective(x))
            end_time = time.perf_counter()

            lst_time.append(end_time - start_time)
            if lst_err != None:
                lst_err[0].append(calc_perc_err(ref, x, ord=1))
                lst_err[1].append(calc_perc_err(ref, x, ord=2))
            save_helper(save, x, k + 1, lst_time, lst_err, idx)

            calc_tol = calc_perc_err(x_old, x, ord=1)
            tolerance.append(calc_tol)
            if verbose:
                pbar.set_postfix(ptol="%0.5f%%" % calc_tol)
                pbar.update()
                pbar.refresh()

            if pAp <= 0 or calc_tol <= ptol:
                break

        if verbose:
            pbar.close()

        if draw_output:
            obj_fig = plt.figure(figsize=(6, 5), dpi=100)
            ax1 = obj_fig.add_subplot(2, 1, 1)
            ax1.plot(np.array(objective), '-', color="red")
            ax1.set_ylabel('objective')
            if lamda is not None:
                ax1.set_title(
                    r'$\frac{1}{2} ||Ax - b||^2_2 + \lambda||x||^2_2 $ using CG')
            else:
                ax1.set_title(
                    r'$\frac{1}{2} ||Ax - b||^2_2 $ using CG')
            ax2 = obj_fig.add_subplot(2, 1, 2, sharex=ax1)
            ax2.plot(np.array(tolerance), '-', color="yellow")
            ax2.set_xlabel('iteration')
            ax2.set_ylabel('tolerance ($x_{i-1}$, $x_i$) (%)')

        return x


def admm(num_iters, ptol, num_normal, A, b, lst_proxg, rho,
         lst_g=None, method="cg", P=None, ref=None,
         save=None, verbose=True, idx=None, draw_output=False):
    r"""ADMM.
    Solves for the following optimization problem:
    .. math::

      \min_x \frac{1}{2} \| A x - y \|_2^2 + \sum_i g_i(x)

    Each iteration involves solving an inner least squares problem which
    is parameterized by the number of A.H * A evaluations.

    Based on:
      For mathematical theory:
      Parikh, N., & Boyd, S.
      Proximal algorithms.
      Foundations and Trends in optimization, 1(3), 127-239.
      DOI: 10.1561/2400000003
      
      And...
      
      For improving convergence times with polynomial preconditioning 
      (a lot of this code borrows Siddarth's impressive polynomial conditioning framework, although this was not disussed in our paper)
      Iyer, S. & Ong, F., et. al.
      Polynomial Preconditioners for Regularized Linear Inverse Problems
      SIAM Journal on Imaging Sciences
      DOI: 10.1137/22M1530355

    Assumes MaxEig(A.H * A) = 1.

    Inputs:
      num_iters (Int): Number of ADMM iterations.
      ptol (Float): l1-percentage tolerance between iterates.
      num_normal (Int): Number of A.H * A evaluations.
      A (Linop): Forward linear operator.
      b (Array): Measurement.
      lst_proxg (List of Prox): List of proximal operators.
      rho (Float): ADMM step size.
      lst_g (List of Functions): List of regularizations. If specified,
                                 objective over iterations are saved.
      method (String): Determines the method used to solve for the inner
                       least squares.
                       - "cg": Conjugate gradient.
                       - "pi": Polynomial inversion.
      P (None or Linop): Preconditioner for "cg".
      save (None or String): If specified, path to save iterations, costs and
                             timings.
      verbose (Bool): Print information.
      idx (Tuple): If passed, slice iterates before saving.
      draw_output (Bool): Plot x after each iteration.
    Returns:
      x (Array): Reconstruction.
    """
    device = sp.get_device(b)
    xp = device.xp

    assert num_normal >= 1
    assert method == "cg" or method == "pi"

    if num_normal == 1 and method == "cg":
        raise Exception("CG requires >= 2 normal evaluations.")

    c = len(lst_proxg)

    def calculate_objective(x):
        obj = sp.to_device(0.5 * xp.linalg.norm(A * x - b)**2, sp.cpu_device)
        for j in range(c):
            obj += sp.to_device(lst_g[j](x), sp.cpu_device)
        return obj

    with device:

        lst_time = []
        lst_err = None
        objective = []
        tolerance = []
        if type(ref) != type(None):
            ref = sp.to_device(ref, device)
            lst_err = ([], [])

        stats = lst_g is not None and save is not None
        if stats:
            lst_obj = []

        start_time = time.perf_counter()
        AHb = A.H(b)
        x = AHb.copy()
        r_AHb = rho * AHb
        xi = xp.zeros([1 + c] + list(x.shape), dtype=xp.complex64)
        ui = xp.zeros_like(xi)
        end_time = time.perf_counter()

        lst_time.append(end_time - start_time)
        if lst_err != None:
            lst_err[0].append(calc_perc_err(ref, x, ord=1))
            lst_err[1].append(calc_perc_err(ref, x, ord=2))
        if stats:
            lst_obj.append(calculate_objective(x))
            save_helper(save, x, 0, lst_time, lst_err, idx, lst_obj)
        else:
            save_helper(save, x, 0, lst_time, lst_err, idx)

        T = sp.linop.Identity(x.shape) + rho * A.N
        if method == "pi":
            P = optpoly.create_polynomial_preconditioner(num_normal - 1, T, 1,
                                                         1 + rho, norm="l_inf",
                                                         verbose=verbose)

        def prox_f(rho, v):
            if method == "cg":
                sp.app.App(sp.alg.ConjugateGradient(A=T, b=v + r_AHb, x=v, P=P,
                                                    max_iter=num_normal - 1,
                                                    tol=1e-4), leave_pbar=False,
                           show_pbar=False,
                           record_time=False).run()
            else:
                v = v - P * (T(v) - (v + r_AHb))
            return v

        if verbose:
            pbar = tqdm(total=num_iters, desc="ADMM", leave=True)

        if draw_output:
            image_fig = plt.figure(figsize=(6, 13.5), dpi=100)

        objective.append(calculate_objective(x))
        for k in range(num_iters):
            x_old = x.copy()
            start_time = time.perf_counter()

            xi[0, ...] = prox_f(rho, x - ui[0, ...])
            for j in range(c):
                xi[1 + j, ...] = lst_proxg[j](rho, x - ui[1 + j, ...])

            x = xp.mean(xi, axis=0)
            ui += xi - x[None, ...]

            if draw_output:
                try:
                    image_slice = x[:, :, np.shape(x)[2]//2]  # 3D
                except:
                    image_slice = x  # 2D
                image_slice_cpu = np.rot90(
                    abs(sp.to_device(image_slice, -1)), 3)

                if num_iters > 144:
                    nrows = 18
                    ncols = 8
                else:
                    nrows = int(np.floor(np.sqrt(num_iters)))
                    ncols = int(np.ceil(np.sqrt(num_iters)))
                if k < 144:
                    ax = image_fig.add_subplot(nrows, ncols, k + 1)
                    ax.axis('off')
                    ax.imshow(image_slice_cpu, cmap="gray")
                    ax.set_title('Iter %d' %
                                 (k + 1), fontsize='x-small', y=0.92)

            objective.append(calculate_objective(x))

            end_time = time.perf_counter()

            lst_time.append(end_time - start_time)
            if lst_err != None:
                lst_err[0].append(calc_perc_err(ref, x, ord=1))
                lst_err[1].append(calc_perc_err(ref, x, ord=2))
            if stats:
                lst_obj.append(calculate_objective(x))
                save_helper(save, x, k + 1, lst_time, lst_err, idx, lst_obj)
            else:
                save_helper(save, x, k + 1, lst_time, lst_err, idx)

            calc_tol = calc_perc_err(x_old, x, ord=1)
            tolerance.append(calc_tol)
            if verbose:
                pbar.set_postfix(ptol="%0.5f%%" % calc_tol)
                pbar.update()
                pbar.refresh()

            if calc_tol <= ptol:
                break

        if draw_output:
            obj_fig = plt.figure(figsize=(6, 5), dpi=100)
            ax1 = obj_fig.add_subplot(2, 1, 1)
            ax1.plot(np.array(objective), '-', color="magenta")
            ax1.set_ylabel('objective')
            if method == "pi":
                ax1.set_title(r'$\frac{1}{2} ||Ax - b||^2_2 + g(x) $ using PI')
            else:
                ax1.set_title(r'$\frac{1}{2} ||Ax - b||^2_2 + g(x) $ using CG')
            ax2 = obj_fig.add_subplot(2, 1, 2, sharex=ax1)
            ax2.plot(np.array(tolerance), '-', color="yellow")
            ax2.set_xlabel('iteration')
            ax2.set_ylabel('tolerance ($x_{i-1}$, $x_i$) (%)')

        if verbose:
            pbar.close()

        return x
