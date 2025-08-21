On the remote machine:

`jupyter lab --no-browser --ip=127.0.0.1 --port 8888`

On the local machine:

`ssh -L 8888:localhost:8888 ncems-no1`

Navigate to `http://localhost:8888` and enter the token from the remote command
