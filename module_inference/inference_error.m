function relerr = inference_error(xhat,x)

relerr = norm(xhat-x,'fro')/norm(x,'fro');

end
        