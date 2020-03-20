# Kinematics Results Database

import os
import sqlite3


class Database:

    def __init__(self, table, memory, columns):
        """
        kinematics propagation database class
        param table: table name
        param memory: use RAM if True else use file for db
        param columns: column names for table
        """
        # NOTE: database file name will be 'table'.db
        # NOTE: dropping db on start to ensure no overlapping results
        self.table = table
        self.db = './results/%s.db' % self.table
        self.drop()
        self.columns = columns
        if memory:
            db_type = ':memory:'
        else:
            db_type = self.db
        self.conn = sqlite3.connect(db_type)
        self.c = self.conn.cursor()
        self.init_db()

    def init_db(self):
        """
        initialise propagation results database
        """
        self.c.execute("""CREATE TABLE IF NOT EXISTS {} ({})"""
                       .format(self.table, self.columns))
        print('...database created\nusing %s table' % self.table)

    def insert(self, *args):
        """
        insert results into propagation database
        """
        n_cols = len(self.columns.split(','))
        n_vars = len(args[0])
        if n_cols != n_vars:
            raise Exception('DB and variable size mismatch')
            return
        mark = '?, ' * n_vars
        self.c.execute("""INSERT INTO {} ({}) VALUES ({})"""
                       .format(self.table, self.columns, mark[:-2]), (args[0]))
        self.conn.commit()

    def query(self, ti=None):
        """
        query table to retrieve results at time
        If time is None, all data from db
        param ti: time to query
        """
        if ti is None:
            command = "SELECT * FROM {}".format(self.table)
            self.c.execute(command)
            print(self.c.fetchall())
        else:
            print('query for time, t = %f' % ti)
            command = "SELECT * FROM {} WHERE t={}".format(self.table, ti)
            self.c.execute(command)
            print(self.columns, '\n', self.c.fetchone())

    def drop(self):
        """
        drop filesystem database
        """
        try:
            os.remove(self.db)
            print('\n... %s database dropped\n' % self.db)
        except OSError:
            print('...database does not exist')

    def close(self):
        """
        close database connection
        """
        print('\n...closing connection\n')
        self.conn.close()
